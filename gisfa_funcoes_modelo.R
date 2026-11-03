################################################################################
# GIS-FA: funcoes auxiliares para salvar/reaplicar modelos finais (sem ASReml)
# A busca MARS continua implementada por mars.forward() do exemplo do usuario;
# aqui apenas reproduzimos as funcoes-base usando os nos exatos salvos.
################################################################################

aplicar_mars_spec <- function(spec, Xpad) {
  Xpad <- as.data.frame(Xpad, check.names = FALSE)
  if (!length(spec)) {
    return(matrix(numeric(0), nrow = nrow(Xpad), ncol = 0L,
                  dimnames = list(rownames(Xpad), NULL)))
  }
  valores <- vector("list", length(spec) + 1L)
  valores[[1L]] <- rep(1, nrow(Xpad))
  for (i in seq_along(spec)) {
    no <- spec[[i]]
    if (!no$variavel %in% names(Xpad)) {
      stop("Falta covariavel MARS para reaplicacao: ", no$variavel)
    }
    if (no$parent_id < 1L || no$parent_id > i) {
      stop("Estrutura MARS invalida: indice do pai no termo ", i)
    }
    z <- Xpad[[no$variavel]]
    hinge <- if (identical(no$direcao, "maior")) {
      pmax(0, z - no$corte)
    } else if (identical(no$direcao, "menor")) {
      pmax(0, no$corte - z)
    } else stop("Direcao MARS invalida")
    valores[[i + 1L]] <- valores[[no$parent_id]] * hinge
  }
  B <- do.call(cbind, valores[-1L])
  dimnames(B) <- list(rownames(Xpad), vapply(spec, `[[`, character(1L), "nome"))
  B
}

criar_bases_mars_final <- function(X, Lambda, mars_forward,
                                   max_features = 10000L,
                                   max_nos = 30L, min_lado = 5L,
                                   min_rss_rel = 0.001,
                                   prop_validacao = 0.20,
                                   seed_validacao = 123L) {
  stopifnot(identical(rownames(X), rownames(Lambda)))
  n <- nrow(X)
  n.valid <- min(max(5L, round(prop_validacao * n)), n - 5L)
  if (n.valid < 1L) stop("Poucos ambientes para dividir FIT/VALID final")
  set.seed(seed_validacao)
  id.valid <- sort(sample(seq_len(n), size = n.valid, replace = FALSE))
  id.fit <- setdiff(seq_len(n), id.valid)
  if (length(id.fit) < 2L * min_lado + 1L) {
    stop("Treinamento final MARS nao possui ambientes FIT suficientes")
  }
  modelos <- list()
  tabelas <- list()
  bases <- list()
  for (f in colnames(Lambda)) {
    medias <- colMeans(X[id.fit, , drop = FALSE])
    desvios <- apply(X[id.fit, , drop = FALSE], 2L, sd)
    keep <- is.finite(desvios) & desvios > 0 & is.finite(medias)
    if (!any(keep)) stop("Nenhuma covariavel MARS valida no fator ", f)
    xpad <- sweep(X[, keep, drop = FALSE], 2L, medias[keep], FUN = "-")
    xpad <- sweep(xpad, 2L, desvios[keep], FUN = "/")
    xpad <- as.data.frame(xpad, check.names = FALSE)
    rownames(xpad) <- rownames(X)
    ajuste <- mars_forward(
      X.train = xpad[id.fit, , drop = FALSE],
      y.train = as.numeric(Lambda[id.fit, f]),
      X.all = xpad,
      max.features = max_features,
      grau.max = 2L,
      min.ambientes.lado = min_lado,
      max.nos = max_nos,
      min.reducao.relativa.RSS = min_rss_rel
    )
    B <- as.matrix(ajuste$B.all)
    if (!ncol(B) || !length(ajuste$spec) || length(ajuste$spec) != ncol(B)) {
      stop("MARS nao produziu funcoes-base portaveis para ", f)
    }
    colunas <- paste0("MARS.", make.names(f), ".BF", seq_len(ncol(B)))
    formulas <- colnames(B)
    colnames(B) <- colunas
    rownames(B) <- rownames(X)
    # Verificacao: o replay das bases deve reproduzir o ajuste numericamente.
    Bagora <- aplicar_mars_spec(ajuste$spec, xpad)
    colnames(Bagora) <- colunas
    if (!isTRUE(all.equal(unname(B), unname(Bagora), tolerance = 1e-10))) {
      stop("Falha na verificacao do replay MARS no fator ", f)
    }
    modelos[[f]] <- list(
      variaveis = colnames(xpad),
      media_fit = medias[keep],
      desvio_fit = desvios[keep],
      spec = ajuste$spec,
      colunas_saida = colunas,
      fator = f, RSS = ajuste$RSS,
      n_fit = length(id.fit), n_valid_interno = length(id.valid)
    )
    bases[[f]] <- B
    tabelas[[f]] <- data.frame(
      fator = f, nome_coluna = colunas, formula = formulas,
      RSS_final = ajuste$RSS, stringsAsFactors = FALSE
    )
  }
  B <- do.call(cbind, bases)
  rownames(B) <- rownames(X)
  desvios <- apply(B, 2L, sd)
  manter <- is.finite(desvios) & desvios > 0
  B <- B[, manter, drop = FALSE]
  if (!ncol(B)) stop("MARS final ficou sem bases com variacao")
  list(X = B, modelos = modelos, colunas_mantidas = colnames(B),
       mapeamento = do.call(rbind, tabelas),
       n_fit = length(id.fit), n_valid_interno = length(id.valid))
}

aplicar_bases_mars_final <- function(X, estrutura) {
  X <- as.matrix(X)
  if (!is.numeric(X) || is.null(colnames(X)) || is.null(rownames(X))) {
    stop("Matriz de covariaveis novas deve ser numerica e nomeada")
  }
  partes <- list()
  for (f in names(estrutura$modelos)) {
    m <- estrutura$modelos[[f]]
    faltam <- setdiff(m$variaveis, colnames(X))
    if (length(faltam)) stop("Covariaveis ausentes no MARS: ", paste(faltam, collapse = ", "))
    xpad <- sweep(X[, m$variaveis, drop = FALSE], 2L,
                  m$media_fit, FUN = "-")
    xpad <- sweep(xpad, 2L, m$desvio_fit, FUN = "/")
    B <- aplicar_mars_spec(m$spec, xpad)
    colnames(B) <- m$colunas_saida
    partes[[f]] <- B
  }
  bases <- do.call(cbind, partes)
  bases <- bases[, estrutura$colunas_mantidas, drop = FALSE]
  if (any(!is.finite(bases))) stop("Bases MARS novas contem NA/Inf")
  bases
}

escolher_configuracoes <- function(results, n_ambientes) {
  tab <- results$summary
  converg <- results$resumo_convergencia
  combos <- unique(tab[, c("factors", "cenario")])
  saida <- list()
  for (i in seq_len(nrow(combos))) {
    k <- combos$factors[i]; cenario <- combos$cenario[i]
    x <- tab[tab$factors == k & tab$cenario == cenario &
               tab$method == "BLUP_fen", , drop = FALSE]
    ok <- converg[converg$FA == k & converg$cenario == cenario &
                    converg$todas_rodadas_concluidas, "ncomp"]
    x <- x[x$ncomp %in% ok & x$n_env.pearson == n_ambientes &
             is.finite(x$cv.pearson), , drop = FALSE]
    if (nrow(x)) {
      x <- x[order(-x$cv.pearson, -x$cv.spearman, x$ncomp,
                   na.last = TRUE), , drop = FALSE]
      saida[[i]] <- data.frame(FA = k, cenario = cenario,
                               ncomp = x$ncomp[1L],
                               cv_pearson = x$cv.pearson[1L],
                               cv_spearman = x$cv.spearman[1L],
                               validacao_completa = TRUE)
    } else {
      saida[[i]] <- data.frame(FA = k, cenario = cenario,
                               ncomp = NA_integer_, cv_pearson = NA_real_,
                               cv_spearman = NA_real_, validacao_completa = FALSE)
    }
  }
  do.call(rbind, saida)
}

salvar_todos_modelos_finais <- function(
    pasta, fa_niveis, cenarios, X, selecao, caracteristicas,
    bienio, periodos, max_comp_pls, mars_forward,
    max_features_mars, max_nos, min_lado, min_rss_rel,
    prop_validacao, seed_validacao) {

  dir.create(pasta, recursive = TRUE, showWarnings = FALSE)
  X <- as.matrix(X)
  if (!is.numeric(X) || is.null(rownames(X)) || is.null(colnames(X)) ||
      any(!is.finite(X))) stop("X de treinamento final invalida")
  status <- list()
  idx <- 0L
  for (k in fa_niveis) {
    fa_path <- file.path(pasta, "FA", paste0("FA", k, "_extraido.rds"))
    if (!file.exists(fa_path)) {
      warning("FA", k, " completo indisponivel, sem modelo final.")
      next
    }
    fa <- readRDS(fa_path)
    nomes <- paste0("fa", seq_len(k))
    Lambda <- as.matrix(fa$rot.loads)[, nomes, drop = FALSE]
    Escores <- as.matrix(fa$rot.scores)[, nomes, drop = FALSE]
    if (!setequal(rownames(Lambda), rownames(X)) ||
        is.null(rownames(Escores)) || any(!is.finite(Lambda)) ||
        any(!is.finite(Escores))) stop("FA final com ambientes ou scores invalidos")
    Lambda <- Lambda[rownames(X), , drop = FALSE]
    estrutura_mars <- NULL
    for (cenario in cenarios) {
      idx <- idx + 1L
      pasta_modelo <- file.path(pasta, cenario, paste0("FA", k))
      dir.create(pasta_modelo, recursive = TRUE, showWarnings = FALSE)
      status_ajuste <- tryCatch({
        if (cenario == "com_MARS") {
          if (is.null(estrutura_mars)) {
            cat("Treinando MARS final FA", k, " com todas as covariaveis...\n")
            estrutura_mars <- criar_bases_mars_final(
              X, Lambda, mars_forward,
              max_features = max_features_mars, max_nos = max_nos,
              min_lado = min_lado, min_rss_rel = min_rss_rel,
              prop_validacao = prop_validacao, seed_validacao = seed_validacao
            )
          }
          Xt <- estrutura_mars$X
          utils::write.table(estrutura_mars$mapeamento,
             file.path(pasta_modelo, "funcoes_MARS_finais.tsv"),
             sep = "\t", quote = FALSE, row.names = FALSE)
        } else {
          Xt <- X
        }
        nmax <- min(nrow(Xt) - 1L, ncol(Xt),
                    qr(scale(Xt))$rank, max_comp_pls)
        if (nmax < 1L) stop("Sem componentes PLS disponiveis")
        modelo_pls <- pls::plsr(
          lambda ~ CovAmb,
          data = data.frame(lambda = I(Lambda), CovAmb = I(Xt)),
          ncomp = nmax, validation = "none", scale = TRUE,
          na.action = "na.fail"
        )
        linha <- selecao[selecao$FA == k & selecao$cenario == cenario, , drop = FALSE]
        selecionado <- if (nrow(linha) && is.finite(linha$ncomp[1L]) &&
                           linha$ncomp[1L] <= nmax) linha$ncomp[1L] else NA_integer_
        if (is.na(selecionado)) {
          warning("FA", k, ", ", cenario,
                  ": sem ncomp validado integralmente. Modelo exportado sem escolha padrao.")
        }
        modelo <- list(
          versao = "GISFA-MARS-portavel-v1", criado_em = Sys.time(),
          fa = k, cenario = cenario,
          ncomp_disponiveis = seq_len(nmax), ncomp_recomendado = selecionado,
          validacao_completa = !is.na(selecionado),
          criterios_LOEO = if (nrow(linha)) linha else NULL,
          PLS = modelo_pls, escores_geneticos = Escores,
          covariaveis_modelo = colnames(X),
          covariaveis_pls = colnames(Xt),
          caracteristicas_diarias = caracteristicas,
          periodos = periodos, bienio_treino = bienio,
          ambiente_treino = rownames(X),
          genotipos = rownames(Escores),
          mars = if (cenario == "com_MARS") estrutura_mars else NULL,
          observacao = paste("Predicao de efeitos geneticos relativos.",
            "Nao corresponde a produtividade absoluta; classificacao relativa", 
            "aos genotipos presentes no treinamento.")
        )
        arquivo <- file.path(pasta_modelo, "modelo_final.rds")
        saveRDS(modelo, arquivo)
        cat("Salvo: ", arquivo, " | ncomp: ",
            if (is.na(selecionado)) "NA" else selecionado, "\n", sep = "")
        data.frame(FA = k, cenario = cenario,
                   arquivo = arquivo, ncomp_max = nmax,
                   ncomp_validado = selecionado, erro = "")
      }, error = function(e) {
        erro <- conditionMessage(e)
        warning("Falhou treino final FA", k, " ", cenario, ": ", erro)
        data.frame(FA = k, cenario = cenario,
                   arquivo = "", ncomp_max = NA_integer_,
                   ncomp_validado = NA_integer_, erro = erro)
      })
      status[[idx]] <- status_ajuste
    }
  }
  if (length(status)) {
    utils::write.table(do.call(rbind, status),
      file.path(pasta, "status_modelos_finais.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  }
  invisible(status)
}

agregar_ciclos_15dias <- function(dados_diarios, bienio, caracteristicas,
                                  periodos_esperados = 1:25) {
  # Entrada: RDS data.frame com lat/lon, YYYYMMDD e colunas ambientais.
  # Agrupamento por local e periodo, media dos 2 ciclos de agosto a julho.
  d <- as.data.frame(dados_diarios)
  essenciais <- c("lat", "lon", "YYYYMMDD", caracteristicas)
  faltam <- setdiff(essenciais, names(d))
  if (length(faltam)) stop("Colunas ausentes no RDS diario: ", paste(faltam, collapse = ", "))
  anos <- suppressWarnings(as.integer(strsplit(bienio, "-", fixed = TRUE)[[1]]))
  if (length(anos) != 2L || anyNA(anos) || anos[2] != anos[1] + 1L) {
    stop("bienio deve ser YYYY-YYYY com anos consecutivos")
  }
  entrada_data <- d$YYYYMMDD
  if (inherits(entrada_data, "Date")) {
    data <- entrada_data
  } else {
    z <- as.character(entrada_data)
    data <- as.Date(z, format = if (all(grepl("^[0-9]{8}$", z[!is.na(z)])))
                          "%Y%m%d" else "%Y-%m-%d")
  }
  if (anyNA(data)) stop("YYYYMMDD tem datas ausentes ou nao reconhecidas")
  ini_ciclo <- as.Date(paste0(ifelse(as.integer(format(data, "%m")) >= 8,
                                       as.integer(format(data, "%Y")),
                                       as.integer(format(data, "%Y")) - 1L), "-08-01"))
  ano_janela <- as.integer(format(data, "%Y")) +
    as.integer(as.integer(format(data, "%m")) >= 8)
  periodo <- as.integer(data - ini_ciclo) %/% 15L + 1L
  d$data_ciclo <- data
  d$ano_janela <- ano_janela
  d$periodo <- periodo
  d <- d[d$ano_janela %in% anos & !is.na(d$lat) & !is.na(d$lon), , drop = FALSE]
  if (!nrow(d)) stop("Nenhuma linha ambiental nos ciclos de ", bienio)
  for (nm in caracteristicas) {
    if (!is.numeric(d[[nm]])) stop("Coluna ambiental nao numerica: ", nm)
  }
  # Mantem o agrupamento original por coordenada; ignora env.x/env.y do RDS.
  chave <- paste(d$lat, d$lon, sep = "|_")
  coords <- unique(data.frame(chave = chave, lat = d$lat, lon = d$lon))
  if (anyDuplicated(coords$chave)) stop("Colisao de chave de coordenadas")
  ids <- sort(unique(chave))
  meta <- coords[match(ids, coords$chave), , drop = FALSE]
  meta$env_novo <- paste0("novo_env.", seq_along(ids))
  # Cobertura: ambos ciclos presentes e cada periodo esperado coberto ao menos uma vez.
  coverage <- expand.grid(chave = ids, ano_janela = anos,
                          periodo = periodos_esperados,
                          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  obs <- unique(data.frame(chave = chave, ano_janela = d$ano_janela,
                           periodo = d$periodo))
  covkey <- paste(coverage$chave, coverage$ano_janela, coverage$periodo)
  obkey <- paste(obs$chave, obs$ano_janela, obs$periodo)
  falta <- coverage[!covkey %in% obkey, , drop = FALSE]
  if (nrow(falta)) {
    detalhes <- utils::head(falta, 12L)
    stop("Cobertura ambiental insuficiente (ambiente, ciclo, periodo): ",
         paste(apply(detalhes, 1L, paste, collapse = ":"), collapse = "; "),
         ". Total combinacoes ausentes: ", nrow(falta))
  }
  grp <- dplyr::group_by(d, lat, lon, periodo)
  avg <- dplyr::summarise(grp,
    dplyr::across(dplyr::all_of(caracteristicas),
      ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)),
    .groups = "drop")
  wide <- tidyr::pivot_wider(avg,
    id_cols = c("lat", "lon"), names_from = "periodo",
    values_from = dplyr::all_of(caracteristicas),
    names_glue = "{.value}.m{periodo}", names_sort = TRUE)
  final <- dplyr::left_join(meta[, c("env_novo", "lat", "lon")], wide,
                            by = c("lat", "lon"))
  final
}

predizer_efeitos_novos <- function(modelo, mat_nova, ncomp = NULL) {
  if (!identical(modelo$versao, "GISFA-MARS-portavel-v1")) {
    stop("Formato de modelo nao reconhecido")
  }
  X <- as.matrix(mat_nova)
  if (!is.numeric(X) || is.null(rownames(X))) stop("Matriz ambiental nova invalida")
  faltam <- setdiff(modelo$covariaveis_modelo, colnames(X))
  if (length(faltam)) stop("Faltam variaveis treinadas: ", paste(faltam, collapse = ", "))
  X <- X[, modelo$covariaveis_modelo, drop = FALSE]
  if (any(!is.finite(X))) stop("Matriz ambiental nova contem NA/Inf nas covariaveis exigidas")
  Z <- if (modelo$cenario == "com_MARS") aplicar_bases_mars_final(X, modelo$mars) else X
  Z <- Z[, modelo$covariaveis_pls, drop = FALSE]
  if (is.null(ncomp)) ncomp <- modelo$ncomp_recomendado
  if (length(ncomp) != 1L || is.na(ncomp) ||
      !ncomp %in% modelo$ncomp_disponiveis) {
    stop("Nao existe ncomp selecionado ou o ncomp solicitado nao esta disponivel")
  }
  preds <- predict(modelo$PLS,
      newdata = data.frame(CovAmb = I(Z)), ncomp = ncomp)
  cargas <- matrix(preds, nrow = nrow(Z), ncol = modelo$fa,
                   dimnames = list(rownames(Z), colnames(modelo$escores_geneticos)))
  efeitos <- cargas %*% t(modelo$escores_geneticos)
  rownames(efeitos) <- rownames(Z)
  list(efeitos = efeitos, cargas = cargas,
       ncomp = ncomp, fa = modelo$fa, cenario = modelo$cenario)
}
