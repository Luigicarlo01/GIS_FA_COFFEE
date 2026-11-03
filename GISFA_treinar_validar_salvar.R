################################################################################
# GIS-FA: comparacao PLS sem MARS vs MARS.forward original + PLS
#
# Fenotipos: media de duas safras por local, cultivar e trial.
# Ambiente: media dos dias do mesmo periodo nos dois ciclos agosto-julho.
# Exemplo: bienio 2024-2025 usa os ciclos encerrados em 2024 e em 2025.
#
# Em CADA rodada LOEO: ambiente excluido ANTES do reajuste FA, do MARS e da PLS.
# FA do treinamento e compartilhada entre os dois cenarios, evitando reajustes
# duplicados. Os fenotipos excluidos sao usados apenas na avaliacao.
# O FA completo e ajustado somente DEPOIS de terminar todas as predicoes LOEO
# para BLUP_BLUP (referencia, nao validacao externa independente).
################################################################################

rm(list = ls())

library(dplyr)
library(asreml)
library(ggplot2)
library(tidyr)
library(tibble)
library(lubridate)
library(knitr)
library(pls)

#==============================================================================
# CONFIGURACAO PBS / ORGANIZACAO DAS SAIDAS
#==============================================================================
WORK_DIR <- Sys.getenv("PBS_O_WORKDIR", unset = getwd())
RESULT_DIR <- Sys.getenv("RESULT_DIR", unset = file.path(WORK_DIR, "resultados_GISFA_15dias_MARS"))
INPUT_DIR <- file.path(RESULT_DIR, "inputs")
SCRIPTS_DIR <- file.path(RESULT_DIR, "scripts_baixados")
LOG_DIR <- file.path(RESULT_DIR, "logs")
DIAG_DIR <- file.path(RESULT_DIR, "diagnosticos")
CHECKPOINT_DIR <- file.path(RESULT_DIR, "checkpoints")
COMP_DIR <- file.path(RESULT_DIR, "comparacao")
CENARIOS <- c("sem_MARS", "com_MARS")
FA_NIVEIS <- c(2L, 3L)

for (pasta in c(RESULT_DIR, INPUT_DIR, SCRIPTS_DIR, LOG_DIR,
                DIAG_DIR, CHECKPOINT_DIR, COMP_DIR,
                file.path(RESULT_DIR, "com_MARS", "variaveis"),
                file.path(RESULT_DIR, "modelos_finais"))) {
  dir.create(pasta, recursive = TRUE, showWarnings = FALSE)
}
for (cenario in CENARIOS) {
  for (k in FA_NIVEIS) {
    dir.create(file.path(RESULT_DIR, cenario, paste0("FA", k)),
               recursive = TRUE, showWarnings = FALSE)
  }
}
setwd(WORK_DIR)
cat("WORK_DIR:", WORK_DIR, "\n")
cat("RESULT_DIR:", RESULT_DIR, "\n")

# Busca o fenotipo no diretorio de submissao (ou pasta inputs).
# Nao substituir silenciosamente data_clean.txt pelo dataset antigo data_ciclo.xlsx.
arquivo.data <- c(
  file.path(WORK_DIR, "data_clean.txt"),
  file.path(WORK_DIR, "inputs", "data_clean.txt"),
  file.path(INPUT_DIR, "data_clean.txt")
)
arquivo.data <- arquivo.data[file.exists(arquivo.data)]
if (!length(arquivo.data)) {
  stop("Falta data_clean.txt. Coloque-o no diretorio de submissao PBS ou em ",
       file.path(WORK_DIR, "inputs", "data_clean.txt"), call. = FALSE)
}
arquivo.data <- normalizePath(arquivo.data[[1L]])
if (normalizePath(dirname(arquivo.data)) != normalizePath(INPUT_DIR)) {
  file.copy(arquivo.data, file.path(INPUT_DIR, "data_clean.txt"), overwrite = TRUE)
}

# Baixa somente se nao houver uma copia local. Se o servidor nao tiver
# acesso ao GitHub, copie antecipadamente o RDS e fa_outs_mod.R para inputs/.
obter_arquivo <- function(nome, url, destino = INPUT_DIR) {
  diretorios <- unique(c(file.path(WORK_DIR, "inputs"), WORK_DIR, destino))
  existentes <- file.path(diretorios, nome)
  existentes <- existentes[file.exists(existentes)]
  alvo <- file.path(destino, nome)
  if (length(existentes)) {
    fonte <- normalizePath(existentes[[1L]])
    if (fonte != normalizePath(alvo, mustWork = FALSE)) {
      if (!file.copy(fonte, alvo, overwrite = TRUE)) {
        stop("Nao foi possivel copiar ", nome, " para ", alvo)
      }
    }
  } else {
    cat("Baixando ", nome, "...\n", sep = "")
    tryCatch(utils::download.file(url, destfile = alvo, mode = "wb", quiet = TRUE),
             error = function(e) stop(
               "Falha ao baixar ", nome, ". Copie este arquivo para ",
               file.path(WORK_DIR, "inputs"), ". Detalhe: ", conditionMessage(e),
               call. = FALSE
             ))
  }
  if (!file.exists(alvo) || file.info(alvo)$size <= 0) {
    stop("Arquivo de entrada vazio/ausente: ", alvo)
  }
  alvo
}

#==============================================================================
# 1. Carregar e organizar os dados fenotipicos
#==============================================================================

data_clean <- read.delim(arquivo.data)
data <- subset(data_clean, yield > 0)

data2 <- data %>%
  group_by(city, lat, lon, cultivar, season, trial) %>%
  summarise(prod = mean(yield, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    env = paste0("env.", as.numeric(as.factor(paste(lat, lon, sep = "_")))),
    env = as.factor(env),
    gen = as.factor(cultivar),
    season = as.factor(season)
  )

data2 <- data2 %>%
  group_by(env) %>%
  filter(n_distinct(season, na.rm = TRUE) > 1) %>%
  ungroup()

resumo <- data2 %>%
  group_by(env) %>%
  summarise(n_anos = n_distinct(season, na.rm = TRUE), .groups = "drop")
print(knitr::kable(resumo, col.names = c("Ambiente", "Numero de anos")))

data3 <- data2 %>%
  mutate(season = as.integer(as.character(season))) %>%
  group_by(env) %>%
  filter(
    n_distinct(season, na.rm = TRUE) != 3L |
      season %in% tail(sort(unique(season[!is.na(season)])), 2L)
  ) %>%
  ungroup() %>%
  mutate(anoN = season)

resumo2 <- data3 %>%
  group_by(env) %>%
  summarise(n_anos = n_distinct(season, na.rm = TRUE), .groups = "drop")
print(knitr::kable(resumo2, col.names = c("Ambiente", "Numero de anos")))

#==============================================================================
# 2. Formar os bienios fenotipicos (inalterado)
#==============================================================================

data3_2021 <- data3 %>%
  mutate(ano_bienio = as.integer(season)) %>%
  filter(ano_bienio >= 2021L)

df_bienio <- data3_2021 %>%
  mutate(par = (ano_bienio - 2022L) %/% 2L) %>%
  group_by(env, par, cultivar, trial) %>%
  filter(n_distinct(ano_bienio) == 2L) %>%
  summarise(
    bienio = paste(sort(unique(ano_bienio)), collapse = "-"),
    across(
      where(is.numeric) &
        !any_of(c("ano", "anoN", "season", "ano_bienio")),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  select(-par) %>%
  mutate(
    env = as.factor(paste(trial, env, sep = "-")),
    cultivar = as.factor(cultivar),
    gen = cultivar,
    bienio = as.factor(bienio)
  )

print(knitr::kable(
  df_bienio %>%
    group_by(trial) %>%
    summarise(
      bienios = paste(sort(unique(as.character(bienio))), collapse = ", "),
      .groups = "drop"
    ),
  col.names = c("Trial", "Bienios disponiveis")
))

#==============================================================================
# 3. Selecionar o bienio e cultivares comuns (inalterado)
#==============================================================================

df_bienio_filter <- subset(df_bienio, bienio == "2024-2025")
df_bienio_filter_sul <- subset(df_bienio_filter, trial == "sul")
df_bienio_filter_state <- subset(df_bienio_filter, trial == "state")

cat("Ambientes sul:", n_distinct(df_bienio_filter_sul$env), "\n")
cat("Ambientes state:", n_distinct(df_bienio_filter_state$env), "\n")

cultivares_comuns <- sort(intersect(
  as.character(df_bienio_filter_sul$cultivar),
  as.character(df_bienio_filter_state$cultivar)
), na.last = NA)
cat("Cultivares comuns:", length(cultivares_comuns), "\n")

df_filter_cultivars <- df_bienio_filter %>%
  filter(
    trial %in% c("state", "sul"),
    cultivar %in% cultivares_comuns,
    !env %in% c("sul-env.42", "sul-env.50")
  ) %>%
  as.data.frame() %>%
  droplevels()

#==============================================================================
# 4. Carregar dados ambientais (inalterado)
#==============================================================================

arquivo.env <- obter_arquivo(
  "new_data_tested_env_diario.RDS",
  paste0("https://raw.githubusercontent.com/Luigicarlo01/GIS_FA_COFFEE/",
         "main/new_data_tested_env_diario.RDS")
)
env.info <- readRDS(arquivo.env) %>%
  select(-any_of(c(
    "FROST_DAYS", "FRUE", "env", "env.x", "env.y", "data", "ano_janela"
  )))

caracteristicas <- env.info %>%
  select(-any_of(c("lat", "lon", "DOY", "YYYYMMDD"))) %>%
  select(where(is.numeric)) %>%
  names()

ambientes <- df_filter_cultivars %>%
  transmute(
    env = as.character(env), lat, lon,
    bienio = as.character(bienio)
  ) %>%
  distinct()

#==============================================================================
# 5. Janelas de 15 dias (inalterado)
#==============================================================================

env.info <- env.info %>%
  mutate(
    data = as.Date(YYYYMMDD),
    ano_janela = ifelse(month(data) >= 8, year(data) + 1, year(data)),
    inicio_ciclo = as.Date(paste0(ano_janela - 1, "-08-01")),
    dia_ciclo = as.integer(data - inicio_ciclo) + 1L,
    periodo = (dia_ciclo - 1L) %/% 15L + 1L
  )

env.list <- list()
for (env_i in unique(ambientes$env)) {
  dt <- ambientes %>% filter(env == env_i)
  anos <- sort(unique(as.integer(unlist(
    strsplit(dt$bienio, "-", fixed = TRUE)
  ))))
  env.frac <- env.info %>%
    semi_join(dt, by = c("lat", "lon")) %>%
    filter(!is.na(data), ano_janela %in% anos) %>%
    group_by(lat, lon, periodo) %>%
    summarise(
      across(
        all_of(caracteristicas),
        ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    mutate(env = env_i) %>%
    relocate(env, lat, lon, periodo)
  env.list[[env_i]] <- env.frac
}

env.info_15dias <- bind_rows(env.list) %>% arrange(env, periodo)
env.final <- env.info_15dias %>%
  pivot_wider(
    id_cols = c(env, lat, lon), names_from = periodo,
    values_from = all_of(caracteristicas),
    names_glue = "{.value}.m{periodo}", names_sort = TRUE
  ) %>%
  arrange(env)

#==============================================================================
# 6. Covariaveis e fenotipos nos mesmos ambientes (inalterado)
#==============================================================================

if (anyDuplicated(env.final$env)) {
  stop("Ha mais de uma linha ambiental para o mesmo env.")
}
matcovamb <- env.final %>%
  dplyr::select(-env, -lat, -lon) %>%
  as.data.frame()
rownames(matcovamb) <- as.character(env.final$env)
CovAmb0 <- as.matrix(matcovamb)
if (!is.numeric(CovAmb0) || !ncol(CovAmb0)) {
  stop("As covariaveis ambientais precisam ser numericas.")
}

ambientes_fenotipicos <- unique(as.character(df_filter_cultivars$env))
ambientes_sem_covariaveis <- setdiff(ambientes_fenotipicos, rownames(CovAmb0))
if (length(ambientes_sem_covariaveis)) {
  warning(
    "Ambientes sem dados ambientais para este bienio: ",
    paste(ambientes_sem_covariaveis, collapse = ", "), call. = FALSE
  )
}
envs <- intersect(ambientes_fenotipicos, rownames(CovAmb0))
if (length(envs) < 3L) {
  stop("Sao necessarios pelo menos tres ambientes com dados ambientais.")
}
CovAmb0 <- CovAmb0[envs, , drop = FALSE]

ok <- apply(CovAmb0, 2L, function(x) {
  if (!all(is.finite(x))) return(FALSE)
  desvio <- stats::sd(x)
  is.finite(desvio) && desvio > 0
})
covariaveis_removidas <- colnames(CovAmb0)[!ok]
if (!any(ok)) stop("Nao restaram covariaveis completas com variacao.")
CovAmb0 <- CovAmb0[, ok, drop = FALSE]

dados_cv <- df_filter_cultivars %>%
  dplyr::filter(as.character(env) %in% envs) %>%
  dplyr::select(env, gen, prod) %>%
  as.data.frame() %>%
  droplevels()
CovAmb0 <- CovAmb0[levels(dados_cv$env), , drop = FALSE]
matcovamb <- as.data.frame(CovAmb0)
stopifnot(
  identical(levels(dados_cv$env), rownames(CovAmb0)),
  all(is.finite(dados_cv$prod)),
  all(is.finite(CovAmb0))
)
cat("Ambientes na validacao:", nrow(CovAmb0), "\n")
cat("Colunas ambientais removidas:", length(covariaveis_removidas), "\n")
cat("Covariaveis usadas antes do MARS:", ncol(CovAmb0), "\n")

#==============================================================================
# 7. LOEO: PLS original e MARS.forward original + PLS, sem vazamento
#==============================================================================

arquivo.fa.outs <- obter_arquivo(
  "fa_outs_mod.R",
  paste0("https://raw.githubusercontent.com/Guilherme-Bravim/possibly_useful/",
         "refs/heads/main/fa_outs_mod.R"),
  destino = SCRIPTS_DIR
)
source(arquivo.fa.outs)

# FUNCOES MARS EXTRAIDAS LITERALMENTE DO EXEMPLO FORNECIDO
# Abaixo preservamos a busca forward, os nos candidatos, os pares de
# funcoes hinge, o RSS e as interacoes de grau ate 2 do exemplo.
# A unica mudanca metodologica e o encaixe dentro de cada treino LOEO.

imputa.media <- function(x) {
  valido <- is.finite(x)

  if (!any(valido)) {
    return(rep(NA_real_, length(x)))
  }

  media <- mean(x[valido], na.rm = TRUE)
  x[!valido] <- media
  x
}

gera.nos.candidatos <- function(
    x,
    suporte,
    min.ambientes.lado = 5,
    max.nos = 30
) {

  x.suporte <- x[suporte & is.finite(x)]
  valores <- sort(unique(x.suporte))

  if (length(valores) < 2 * min.ambientes.lado + 1) {
    return(numeric(0))
  }

  candidatos <- valores[
    vapply(
      valores,
      function(k) {
        sum(x.suporte < k) >= min.ambientes.lado &&
          sum(x.suporte > k) >= min.ambientes.lado
      },
      logical(1)
    )
  ]

  if (length(candidatos) == 0) {
    return(numeric(0))
  }

  if (is.finite(max.nos) && length(candidatos) > max.nos) {
    posicoes <- unique(
      round(
        seq(
          1,
          length(candidatos),
          length.out = max.nos
        )
      )
    )
    candidatos <- candidatos[posicoes]
  }

  candidatos
}

mars.forward <- function(
    X.train,
    y.train,
    X.all,
    max.features = 40,
    grau.max = 2,
    min.ambientes.lado = 5,
    max.nos = 30,
    min.reducao.relativa.RSS = 0.001
) {

  X.train <- as.data.frame(X.train)
  X.all <- as.data.frame(X.all)
  y.train <- as.numeric(y.train)

  n.train <- nrow(X.train)
  n.all <- nrow(X.all)

  B.train <- matrix(1, nrow = n.train, ncol = 1)
  B.all <- matrix(1, nrow = n.all, ncol = 1)

  colnames(B.train) <- "(Intercept)"
  colnames(B.all) <- "(Intercept)"

  termos <- list(
    list(
      nome = "(Intercept)",
      grau = 0L,
      variaveis = character(0),
      base.train = rep(1, n.train),
      base.all = rep(1, n.all),
      parent_id = NA_integer_, variavel = NA_character_,
      corte = NA_real_, direcao = NA_character_
    )
  )

  ajuste.atual <- lm.fit(x = B.train, y = y.train)
  rss.atual <- sum(ajuste.atual$residuals^2)

  max.pares.features <- floor(max.features / 2)
  max.pares.amostras <- floor((n.train - 1) / 2)
  max.pares <- min(max.pares.features, max.pares.amostras)

  if (
    max.pares < 1 ||
      !is.finite(rss.atual) ||
      rss.atual <= 0
  ) {
    return(
      list(
        B.all = X.all[, FALSE, drop = FALSE],
        RSS = rss.atual,
        spec = list()
      )
    )
  }

  nomes.selecionados <- "(Intercept)"

  for (passo in seq_len(max.pares)) {

    if (
      !is.finite(rss.atual) ||
        rss.atual <= .Machine$double.eps
    ) {
      break
    }

    melhor <- NULL
    melhor.rss <- rss.atual

    for (indice.pai in seq_along(termos)) {

      pai <- termos[[indice.pai]]

      if (pai$grau >= grau.max) {
        next
      }

      suporte <- abs(pai$base.train) > sqrt(.Machine$double.eps)

      if (sum(suporte) < 2 * min.ambientes.lado + 1) {
        next
      }

      for (variavel in names(X.train)) {

        if (variavel %in% pai$variaveis) {
          next
        }

        x.train <- X.train[[variavel]]
        x.all <- X.all[[variavel]]

        nos <- gera.nos.candidatos(
          x = x.train,
          suporte = suporte,
          min.ambientes.lado = min.ambientes.lado,
          max.nos = max.nos
        )

        if (length(nos) == 0) {
          next
        }

        for (k in nos) {

          base.maior.train <- pai$base.train * pmax(0, x.train - k)
          base.menor.train <- pai$base.train * pmax(0, k - x.train)

          desvio.maior <- sd(base.maior.train)
          desvio.menor <- sd(base.menor.train)

          if (
            !is.finite(desvio.maior) ||
              !is.finite(desvio.menor) ||
              desvio.maior <= 0 ||
              desvio.menor <= 0
          ) {
            next
          }

          k.texto <- format(signif(k, 8), trim = TRUE)

          if (pai$nome == "(Intercept)") {
            prefixo <- ""
          } else {
            prefixo <- paste0("(", pai$nome, ")*")
          }

          nome.maior <- paste0(
            prefixo,
            "h(",
            variavel,
            "-",
            k.texto,
            ")"
          )

          nome.menor <- paste0(
            prefixo,
            "h(",
            k.texto,
            "-",
            variavel,
            ")"
          )

          if (
            nome.maior %in% nomes.selecionados ||
              nome.menor %in% nomes.selecionados
          ) {
            next
          }

          B.candidato <- cbind(
            B.train,
            base.maior.train,
            base.menor.train
          )

          ajuste.candidato <- lm.fit(
            x = B.candidato,
            y = y.train
          )

          if (ajuste.candidato$rank < ncol(B.candidato)) {
            next
          }

          rss.candidato <- sum(ajuste.candidato$residuals^2)

          if (
            is.finite(rss.candidato) &&
              rss.candidato < melhor.rss
          ) {

            melhor.rss <- rss.candidato

            melhor <- list(
              parent_id = indice.pai,
              variavel = variavel, corte = k,
              nome.maior = nome.maior,
              nome.menor = nome.menor,
              grau = pai$grau + 1L,
              variaveis = c(
                pai$variaveis,
                variavel
              ),
              base.maior.train = base.maior.train,
              base.menor.train = base.menor.train,
              base.maior.all = pai$base.all * pmax(
                0,
                x.all - k
              ),
              base.menor.all = pai$base.all * pmax(
                0,
                k - x.all
              ),
              RSS = rss.candidato
            )
          }
        }
      }
    }

    if (is.null(melhor)) {
      break
    }

    reducao.relativa.RSS <- (
      rss.atual - melhor$RSS
    ) / rss.atual

    if (reducao.relativa.RSS < min.reducao.relativa.RSS) {
      break
    }

    B.train <- cbind(
      B.train,
      melhor$base.maior.train,
      melhor$base.menor.train
    )

    colnames(B.train)[
      (ncol(B.train) - 1):ncol(B.train)
    ] <- c(
      melhor$nome.maior,
      melhor$nome.menor
    )

    B.all <- cbind(
      B.all,
      melhor$base.maior.all,
      melhor$base.menor.all
    )

    colnames(B.all)[
      (ncol(B.all) - 1):ncol(B.all)
    ] <- c(
      melhor$nome.maior,
      melhor$nome.menor
    )

    termos[[length(termos) + 1]] <- list(
      nome = melhor$nome.maior,
      grau = melhor$grau,
      variaveis = melhor$variaveis,
      base.train = melhor$base.maior.train,
      base.all = melhor$base.maior.all,
      parent_id = melhor$parent_id, variavel = melhor$variavel,
      corte = melhor$corte, direcao = "maior"
    )

    termos[[length(termos) + 1]] <- list(
      nome = melhor$nome.menor,
      grau = melhor$grau,
      variaveis = melhor$variaveis,
      base.train = melhor$base.menor.train,
      base.all = melhor$base.menor.all,
      parent_id = melhor$parent_id, variavel = melhor$variavel,
      corte = melhor$corte, direcao = "menor"
    )

    nomes.selecionados <- c(
      nomes.selecionados,
      melhor$nome.maior,
      melhor$nome.menor
    )

    rss.atual <- melhor$RSS
  }

  B.final <- B.all[
    ,
    -1,
    drop = FALSE
  ]

  B.final <- as.data.frame(
    B.final,
    check.names = FALSE
  )

  # Metadados adicionais, SEM modificar a busca forward do exemplo.
  # Guarda os nos sem arredondamento para reproduzir as bases em ambientes novos.
  specs <- lapply(termos[-1L], function(t) list(
    nome = t$nome, parent_id = t$parent_id,
    variavel = t$variavel, corte = t$corte, direcao = t$direcao
  ))
  list(
    B.all = B.final,
    RSS = rss.atual,
    spec = specs
  )
}

cross.val.loeo <- function(
    factors = c(2L, 3L), phen.data, env.data, n.comp.pls = 2:10,
    description = "LOEO_retreino_MARS", out.dir = ".", fa.outs.fun = fa.outs,
    overwrite = FALSE, verbose = TRUE, maxit = 1000L, max.updates = 2L,
    compare.full.blups = TRUE,
    # Parametros ORIGINAIS do mars.forward() do exemplo enviado.
    cenarios = c("sem_MARS", "com_MARS"),
    max.features.por.FA = 10000L,
    max.nos.candidatos = 30L,
    min.ambientes.lado = 5L,
    min.reducao.relativa.RSS = 0.001,
    usar.interacoes = TRUE,
    prop.validacao = 0.20,
    seed.validacao = 123L
) {
  cenarios <- unique(cenarios)
  if (!length(cenarios) || !all(cenarios %in% c("sem_MARS", "com_MARS"))) {
    stop("cenarios deve conter 'sem_MARS' e/ou 'com_MARS'.")
  }
  for (pkg in c("asreml", "pls", "dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Instale o pacote: ", pkg)
  }
  check_int <- function(x, nome, minimo = 1L) {
    if (!is.numeric(x) || !length(x) || anyNA(x) ||
        any(!is.finite(x) | x < minimo | x != floor(x))) {
      stop(nome, " deve conter inteiros >= ", minimo, ".")
    }
    sort(unique(as.integer(x)))
  }
  factors <- check_int(factors, "factors")
  n.comp.pls <- check_int(n.comp.pls, "n.comp.pls")
  maxit <- check_int(maxit, "maxit")
  max.updates <- check_int(max.updates, "max.updates", minimo = 0L)
  max.features.por.FA <- check_int(max.features.por.FA, "max.features.por.FA")
  max.nos.candidatos <- check_int(max.nos.candidatos, "max.nos.candidatos")
  min.ambientes.lado <- check_int(min.ambientes.lado, "min.ambientes.lado")
  seed.validacao <- check_int(seed.validacao, "seed.validacao", minimo = 0L)
  if (!is.logical(usar.interacoes) || length(usar.interacoes) != 1L ||
      is.na(usar.interacoes)) stop("usar.interacoes deve ser TRUE/FALSE.")
  if (length(maxit) != 1L || length(max.updates) != 1L ||
      length(max.features.por.FA) != 1L || length(max.nos.candidatos) != 1L ||
      length(min.ambientes.lado) != 1L || length(seed.validacao) != 1L) {
    stop("Os parametros de controle devem ser escalares.")
  }
  if (!is.numeric(prop.validacao) || length(prop.validacao) != 1L ||
      !is.finite(prop.validacao) || prop.validacao <= 0 || prop.validacao >= 1) {
    stop("prop.validacao precisa estar entre 0 e 1.")
  }
  if (!is.numeric(min.reducao.relativa.RSS) ||
      length(min.reducao.relativa.RSS) != 1L ||
      !is.finite(min.reducao.relativa.RSS) ||
      min.reducao.relativa.RSS < 0) {
    stop("min.reducao.relativa.RSS precisa ser >= 0.")
  }
  grau.max <- if (usar.interacoes) 2L else 1L
  if (!is.logical(compare.full.blups) || length(compare.full.blups) != 1L ||
      is.na(compare.full.blups)) {
    stop("compare.full.blups deve ser TRUE ou FALSE.")
  }
  if (!is.function(fa.outs.fun)) stop("fa.outs.fun precisa ser funcao.")
  if (!is.character(description) || length(description) != 1L ||
      is.na(description) || !nzchar(description) ||
      basename(description) != description) {
    stop("description deve ser um nome sem diretorios.")
  }

  phen.data <- as.data.frame(phen.data)
  if (!all(c("env", "gen", "prod") %in% names(phen.data))) {
    stop("phen.data precisa conter env, gen e prod.")
  }
  if (!nrow(phen.data) || anyNA(phen.data$env) || anyNA(phen.data$gen) ||
      !is.numeric(phen.data$prod) || any(!is.finite(phen.data$prod))) {
    stop("phen.data invalido; confira env, gen e prod.")
  }
  phen.data$env <- droplevels(factor(phen.data$env))
  phen.data$gen <- droplevels(factor(phen.data$gen))
  name.env <- levels(phen.data$env)
  name.gen <- levels(phen.data$gen)
  if (length(name.env) < 3L) stop("Sao necessarios pelo menos 3 ambientes.")

  X <- as.matrix(env.data)
  if (!is.numeric(X) || !ncol(X) || is.null(rownames(X)) ||
      anyDuplicated(rownames(X))) {
    stop("env.data precisa ter covariaveis numericas e ambientes unicos nas linhas.")
  }
  if (length(setdiff(name.env, rownames(X)))) {
    stop("Ambientes sem covariaveis: ",
         paste(setdiff(name.env, rownames(X)), collapse = ", "))
  }
  X <- X[name.env, , drop = FALSE]
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  if (anyDuplicated(colnames(X)) || any(!is.finite(X))) {
    stop("Covariaveis duplicadas ou NA/NaN/Inf em env.data.")
  }
  if (any(!nzchar(colnames(X)))) stop("Colunas ambientais sem nome.")

  # Checar a viabilidade da grade PLS em CADA conjunto de treinamento.
  fold.cov <- setNames(lapply(name.env, function(excluido) {
    ids <- setdiff(name.env, excluido)
    tr <- X[ids, , drop = FALSE]
    desvios <- apply(tr, 2, stats::sd)
    keep <- is.finite(desvios) & desvios > 0
    if (!any(keep)) stop("Nao ha covariaveis variaveis sem ", excluido)
    limite <- min(nrow(tr) - 1L, sum(keep),
                  qr(scale(tr[, keep, drop = FALSE]))$rank)
    list(keep = keep, limit = limite)
  }), name.env)
  max_allowed <- min(vapply(fold.cov, `[[`, integer(1), "limit"))
  if (max(n.comp.pls) > max_allowed) {
    stop("A grade PLS requer ate ", max(n.comp.pls),
         " componentes, mas o limite comum e ", max_allowed, ".")
  }

  out.dir <- path.expand(out.dir)
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out.dir)) stop("Nao foi possivel criar out.dir.")
  for (cenario in cenarios) for (k in factors) {
    dir.create(file.path(out.dir, cenario, paste0("FA", k)),
               recursive = TRUE, showWarnings = FALSE)
  }
  for (subdir in c("diagnosticos", "checkpoints", "logs")) {
    dir.create(file.path(out.dir, subdir), recursive = TRUE, showWarnings = FALSE)
  }
  file_base <- function(cenario, k, j) file.path(
    out.dir, cenario, paste0("FA", k),
    paste0(description, "_", cenario, "_FA", k, "_ncomp", j)
  )
  paths <- unlist(lapply(cenarios, function(cenario) unlist(
    lapply(factors, function(k) unlist(lapply(n.comp.pls, function(j)
      paste0(file_base(cenario, k, j),
             c(".txt", "_por_ambiente.txt", "_predicoes.txt")))))
  )))
  if (!overwrite && any(file.exists(paths))) {
    stop("Arquivos ja existentes. Mude description/out.dir ou overwrite = TRUE.")
  }
  write_tab <- function(z, path) utils::write.table(
    z, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
  )
  bind_df <- function(x) if (!length(x)) NULL else do.call(rbind, x)
  safe_cor <- function(a, b, metodo) {
    ok <- is.finite(a) & is.finite(b)
    a <- a[ok]; b <- b[ok]
    if (length(a) < 2L || length(unique(a)) < 2L ||
        length(unique(b)) < 2L) return(NA_real_)
    suppressWarnings(stats::cor(a, b, method = metodo))
  }
  desc_cor <- function(z) {
    z <- z[is.finite(z)]
    n <- length(z)
    s <- if (n > 1L) stats::sd(z) else NA_real_
    c(mean = if (n) mean(z) else NA_real_, sd = s,
      se = if (n > 1L) s / sqrt(n) else NA_real_, n = n)
  }
  pair_metrics <- function(env, metodo, pred, obs, gens) {
    valid <- is.finite(pred) & is.finite(obs)
    data.frame(
      Ambiente = env, method = metodo,
      Pearson_Correlation = safe_cor(pred, obs, "pearson"),
      Spearman_Correlation = safe_cor(pred, obs, "spearman"),
      n_pairs = sum(valid), n_gen = length(unique(gens[valid])),
      stringsAsFactors = FALSE
    )
  }

  # MARS DO EXEMPLO: forward sobre TODAS as covariaveis validas, sem SHAP
  # nem ranking previo. Um ajuste para CADA coluna das cargas FA do treino.
  # Cria pares hinge e interacoes ate grau.max. A PLS usa as FUNCOES-BASE.
  # Divisao interna 80/20: MARS usa os 80% FIT, totalmente dentro do LOEO;
  # PLS usa o treino LOEO completo e o ambiente retido so para predicao.
  construir_bases_mars <- function(Xtreino, Xteste, Lambda, excluido, k) {
    stopifnot(
      identical(rownames(Xtreino), rownames(Lambda)),
      nrow(Xteste) == 1L,
      !excluido %in% rownames(Xtreino)
    )
    n.train <- nrow(Xtreino)
    n.valid <- max(5L, round(prop.validacao * n.train))
    n.valid <- min(n.valid, n.train - 5L)
    if (n.valid < 1L) stop("Nao ha ambientes suficientes para split FIT/VALID.")
    # O mesmo split interno do exemplo, reinicializado por rodada LOEO.
    set.seed(seed.validacao)
    id.valid <- sort(sample(seq_len(n.train), size = n.valid, replace = FALSE))
    id.fit <- setdiff(seq_len(n.train), id.valid)
    if (length(id.fit) < 2L * min.ambientes.lado + 1L) {
      stop("FIT do MARS tem ", length(id.fit), " ambientes; sao necessarios pelo menos ",
           2L * min.ambientes.lado + 1L, " para gerar nos.")
    }

    # Para cada fator, medias e desvios calculados EXCLUSIVAMENTE nos 80% FIT.
    # Os ambientes VALID internos e o excluido usam essa mesma escala.
    Xall <- rbind(Xtreino, Xteste)
    bases_fatores <- list()
    mapas <- list()
    resumos <- list()
    for (f in colnames(Lambda)) {
      medias <- colMeans(Xtreino[id.fit, , drop = FALSE])
      desvios <- apply(Xtreino[id.fit, , drop = FALSE], 2L, stats::sd)
      validas <- is.finite(desvios) & desvios > 0 & is.finite(medias)
      if (!any(validas)) {
        stop("Nenhuma covariavel variavel no FIT MARS: ", f, ", sem ", excluido)
      }
      Xpad <- sweep(Xall[, validas, drop = FALSE], 2L, medias[validas], "-")
      Xpad <- sweep(Xpad, 2L, desvios[validas], "/")
      Xpad <- as.data.frame(Xpad, check.names = FALSE)
      y <- as.numeric(Lambda[, f])
      # A funcao mars.forward() abaixo foi copiada do exemplo original.
      if (verbose) message("    mars.forward: FA", k, ", sem ", excluido,
                           ", fator ", f, ", ", sum(validas), " covariaveis candidatas.")
      inicio_mars <- proc.time()[["elapsed"]]
      ajuste <- mars.forward(
        X.train = Xpad[id.fit, , drop = FALSE],
        y.train = y[id.fit],
        X.all = Xpad,
        max.features = max.features.por.FA,
        grau.max = grau.max,
        min.ambientes.lado = min.ambientes.lado,
        max.nos = max.nos.candidatos,
        min.reducao.relativa.RSS = min.reducao.relativa.RSS
      )
      if (verbose) message("    mars.forward: FA", k, ", sem ", excluido,
                           ", fator ", f, " concluido em ",
                           round(proc.time()[["elapsed"]] - inicio_mars, 1),
                           " s; RSS = ", signif(ajuste$RSS, 6))
      B <- as.matrix(ajuste$B.all)
      if (!ncol(B)) {
        stop("mars.forward nao produziu funcoes-base: ", f, ", sem ", excluido)
      }
      rownames(B) <- rownames(Xall)
      formulas <- colnames(B)
      nomes <- paste0("MARS.", make.names(f), ".BF", seq_len(ncol(B)))
      colnames(B) <- nomes
      bases_fatores[[f]] <- B
      mapas[[f]] <- data.frame(
        factors = k, Ambiente = excluido, fator = f,
        feature_saida = nomes, funcao_MARS = formulas,
        n_candidatas = sum(validas), n_fit = length(id.fit),
        n_valid_interno = length(id.valid),
        stringsAsFactors = FALSE
      )
      resumos[[f]] <- data.frame(
        factors = k, Ambiente = excluido, fator = f,
        n_features_entrada = sum(validas), n_features_MARS = ncol(B),
        RSS_final = ajuste$RSS, n_fit = length(id.fit),
        n_valid_interno = length(id.valid),
        stringsAsFactors = FALSE
      )
    }
    # Como no exemplo, concatenar todas as funcoes-base dos fatores.
    Btotal <- do.call(cbind, bases_fatores)
    colnames(Btotal) <- unlist(lapply(bases_fatores, colnames), use.names = FALSE)
    rownames(Btotal) <- rownames(Xall)
    # Filtrar apenas usando o TREINO LOEO. Nao completar com variaveis
    # originais nem usar fenotipos/cargas do ambiente excluido.
    desvio_treino <- apply(Btotal[rownames(Xtreino), , drop = FALSE], 2L, stats::sd)
    manter <- is.finite(desvio_treino) & desvio_treino > 0
    Btotal <- Btotal[, manter, drop = FALSE]
    if (!ncol(Btotal)) stop("MARS: funcoes-base sem variacao no treino.")
    mapeamento <- do.call(rbind, mapas)
    mapeamento$mantida_para_PLS <- mapeamento$feature_saida %in% colnames(Btotal)
    list(
      train = Btotal[rownames(Xtreino), , drop = FALSE],
      test = Btotal[excluido, , drop = FALSE],
      map = mapeamento, resumo = do.call(rbind, resumos),
      n_fit = length(id.fit), n_valid_interno = length(id.valid)
    )
  }

  fit_attempts <- list()
  fold_diagnostics <- list()
  features <- list()
  mars_resumos <- list()
  summaries <- list()
  by_environment <- list()
  predictions <- list()
  reference_diagnostics <- list()
  pred_current <- list()
  k_atual <- NA_integer_

  # Logs/auditoria: uma linha para cada ambiente, FA, cenario e ncomp.
  convergence_tables <- function() {
    tab <- expand.grid(
      Ambiente_retirado = name.env, FA = factors, cenario = cenarios,
      ncomp = n.comp.pls, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    tab <- tab[order(tab$FA, tab$cenario, tab$ncomp, tab$Ambiente_retirado), ]
    rownames(tab) <- NULL
    tab$FA_convergiu <- NA
    tab$PLS_valida <- FALSE
    tab$MARS_aplicado <- NA
    tab$n_covariaveis <- NA_integer_
    tab$status <- "Nao executado"
    tent <- bind_df(fit_attempts)
    if (!is.null(tent)) {
      ult <- tent[!duplicated(tent$contexto, fromLast = TRUE), ]
      alvo <- paste0("FA", tab$FA, " sem ", tab$Ambiente_retirado)
      idx <- match(alvo, ult$contexto)
      has <- !is.na(idx)
      tab$FA_convergiu[has] <- ult$converged[idx[has]]
      tab$status[has & tab$FA_convergiu %in% FALSE] <- "FA nao convergiu"
      tab$status[has & tab$FA_convergiu %in% TRUE] <- "FA ok; PLS pendente"
    }
    diag <- bind_df(fold_diagnostics)
    if (!is.null(diag)) {
      chave <- paste(tab$Ambiente_retirado, tab$FA, tab$cenario, sep = "\r")
      idx <- match(chave, paste(diag$Ambiente, diag$factors,
                                diag$cenario, sep = "\r"))
      valid <- !is.na(idx)
      tab$PLS_valida[valid] <- diag$PLS_valida[idx[valid]] &
        tab$ncomp[valid] <= diag$ncomp_max[idx[valid]]
      tab$MARS_aplicado[valid] <- diag$MARS_aplicado[idx[valid]]
      tab$n_covariaveis[valid] <- diag$n_covariaveis[idx[valid]]
      tab$status[valid & !tab$PLS_valida] <- "PLS indisponivel neste ncomp"
      tab$status[valid & tab$PLS_valida] <- "Concluido"
    }
    res <- tab %>%
      group_by(FA, cenario, ncomp) %>%
      summarise(
        n_ambientes = n(),
        FA_convergiram = sum(FA_convergiu %in% TRUE),
        PLS_concluidas = sum(PLS_valida),
        todas_rodadas_concluidas = all(status == "Concluido"),
        .groups = "drop"
      ) %>% as.data.frame()
    list(detalhada = tab, resumo = res)
  }
  checkpoint <- function() {
    tb <- convergence_tables()
    write_tab(tb$detalhada,
              file.path(out.dir, "diagnosticos", "convergencia_por_ambiente_FA_cenario_ncomp.txt"))
    write_tab(tb$resumo,
              file.path(out.dir, "diagnosticos", "resumo_convergencia_FA_cenario_ncomp.txt"))
    saveRDS(list(
      summary = bind_df(summaries),
      by_environment = bind_df(by_environment),
      predictions = predictions,
      diagnostics = bind_df(fold_diagnostics),
      mars_features = bind_df(features),
      mars_resumo = bind_df(mars_resumos),
      fit_attempts = bind_df(fit_attempts),
      reference_diagnostics = bind_df(reference_diagnostics),
      convergencia = tb$detalhada,
      resumo_convergencia = tb$resumo,
      current_factors = k_atual,
      current_predictions = pred_current
    ), file.path(out.dir, "checkpoints", "checkpoint_LOEO_PLS_MARS.rds"))
  }
  on.exit(tryCatch(checkpoint(), error = function(e) warning(
    "Falha no checkpoint: ", conditionMessage(e), call. = FALSE
  )), add = TRUE)

  for (k in factors) {
    k_atual <- k
    if (verbose) message("FA", k, ": iniciando rodadas LOEO.")
    random_formula <- as.formula(paste0("~ gen:fa(env, ", k, ")"))

    fit_fa <- function(dat, contexto, tipo_ajuste = "LOEO") {
      dat <- droplevels(as.data.frame(dat))
      dat <- dat[order(dat$env, dat$gen), , drop = FALSE]
      controles <- data.frame(
        step.size = c(0.1, 0.316, 0.05),
        ai.penalty = c(20, 10, 30)
      )
      mod <- NULL
      ultimo_mod <- NULL
      ultimo_erro <- "O ajuste terminou sem convergencia."
      for (attempt in seq_len(max.updates + 1L)) {
        perfil <- (attempt - 1L) %% nrow(controles) + 1L
        passo <- controles$step.size[perfil]
        penalizacao <- controles$ai.penalty[perfil]
        reiniciar <- is.null(mod) || perfil == 3L
        avisos <- character(0)
        inicio <- proc.time()[["elapsed"]]
        cand <- tryCatch(withCallingHandlers({
          if (reiniciar) {
            asreml::asreml(
              fixed = prod ~ env, random = random_formula,
              data = dat, maxit = maxit,
              step.size = passo, ai.penalty = penalizacao,
              na.action = asreml::na.method(x = "include", y = "include")
            )
          } else {
            stats::update(mod, maxit = maxit,
                          step.size = passo, ai.penalty = penalizacao)
          }
        }, warning = function(w) {
          avisos <<- c(avisos, conditionMessage(w))
          invokeRestart("muffleWarning")
        }), error = function(e) e)
        erro <- inherits(cand, "error")
        convergiu <- !erro && isTRUE(cand$converge)
        loglik <- if (!erro && is.numeric(cand$loglik) &&
                      length(cand$loglik) == 1L) cand$loglik else NA_real_
        fit_attempts[[length(fit_attempts) + 1L]] <<- data.frame(
          factors = k, contexto = contexto, tipo_ajuste = tipo_ajuste,
          tentativa = attempt,
          tempo_segundos = proc.time()[["elapsed"]] - inicio,
          inicializacao = if (reiniciar) "nova" else "continuacao",
          step.size = passo, ai.penalty = penalizacao,
          n_env = nlevels(dat$env), n_gen = nlevels(dat$gen), n_obs = nrow(dat),
          converged = convergiu, loglik = loglik,
          avisos = gsub("[\r\n\t]+", " ", paste(unique(avisos), collapse = " | ")),
          erro = if (erro) gsub("[\r\n\t]+", " ", conditionMessage(cand)) else "",
          stringsAsFactors = FALSE
        )
        write_tab(bind_df(fit_attempts),
                  file.path(out.dir, "logs", "tentativas_ajuste_FA.txt"))
        if (verbose) message(contexto, ": tentativa ", attempt,
                             "; convergiu = ", convergiu,
                             "; logLik = ", signif(loglik, 7))
        if (convergiu) return(cand)
        if (erro) {
          ultimo_erro <- conditionMessage(cand)
          mod <- NULL
        } else {
          mod <- cand
          ultimo_mod <- cand
          ultimo_erro <- "O ajuste terminou sem convergencia."
        }
      }
      if (!is.null(ultimo_mod)) {
        prefix <- file.path(out.dir, "diagnosticos", paste0(
          gsub("[^[:alnum:]_.-]+", "_", contexto), "_nao_convergido"
        ))
        tryCatch({
          saveRDS(ultimo_mod, paste0(prefix, ".rds"))
          vc <- as.data.frame(summary(ultimo_mod)$varcomp)
          vc <- data.frame(componente = rownames(vc), vc,
                           row.names = NULL, check.names = FALSE)
          write_tab(vc, paste0(prefix, "_componentes.txt"))
        }, error = function(e) warning(conditionMessage(e), call. = FALSE))
      }
      stop(contexto, ": ", ultimo_erro,
           ". Consulte tentativas_ajuste_FA.txt e arquivos *_nao_convergido.",
           call. = FALSE)
    }

    # Matrizes separadas por cenario/ncomp, SEM dividir o FA do treino.
    pred_current <- setNames(lapply(cenarios, function(cenario) {
      lapply(n.comp.pls, function(j) matrix(
        NA_real_, nrow = length(name.gen), ncol = length(name.env),
        dimnames = list(name.gen, name.env)
      ))
    }), cenarios)

    for (excluido in name.env) {
      if (verbose) message("FA", k, ": sem ", excluido)
      dat.train <- droplevels(
        phen.data[as.character(phen.data$env) != excluido, , drop = FALSE]
      )
      dat.train <- dat.train[order(dat.train$env), , drop = FALSE]
      mod <- fit_fa(dat.train, paste0("FA", k, " sem ", excluido))
      fa <- fa.outs.fun(mod, name.env = "env", name.gen = "gen")
      nomes_fa <- paste0("fa", seq_len(k))
      Lambda <- as.matrix(fa$rot.loads)
      Escores <- as.matrix(fa$rot.scores)
      if (!all(nomes_fa %in% colnames(Lambda)) ||
          !all(nomes_fa %in% colnames(Escores))) {
        stop("FA", k, " sem ", excluido,
             ": cargas/escores sem os fatores esperados.")
      }
      Lambda <- Lambda[, nomes_fa, drop = FALSE]
      Escores <- Escores[, nomes_fa, drop = FALSE]
      if (is.null(rownames(Lambda)) || anyDuplicated(rownames(Lambda)) ||
          !setequal(rownames(Lambda), levels(dat.train$env)) ||
          is.null(rownames(Escores)) || anyDuplicated(rownames(Escores)) ||
          any(!is.finite(Lambda)) || any(!is.finite(Escores))) {
        stop("FA", k, " sem ", excluido, ": cargas/escores invalidos.")
      }
      keep <- fold.cov[[excluido]]$keep
      Xtrain <- X[rownames(Lambda), keep, drop = FALSE]
      Xtest <- X[excluido, keep, drop = FALSE]
      stopifnot(
        !excluido %in% rownames(Lambda),
        identical(rownames(Lambda), rownames(Xtrain))
      )

      bases_mars <- NULL
      erro_mars <- ""
      if ("com_MARS" %in% cenarios) {
        bases_mars <- tryCatch(
          construir_bases_mars(Xtrain, Xtest, Lambda, excluido, k),
          error = function(e) {
            erro_mars <<- conditionMessage(e)
            warning("MARS falhou no FA", k, " sem ", excluido,
                    ": ", erro_mars, call. = FALSE)
            NULL
          }
        )
        if (!is.null(bases_mars)) {
          features[[paste(k, excluido, sep = "|")]] <- bases_mars$map
          mars_resumos[[paste(k, excluido, sep = "|")]] <- bases_mars$resumo
          if (verbose) message(
            "  MARS original: ", ncol(bases_mars$train),
            " funcoes-base para PLS; ", bases_mars$n_fit, " ambientes FIT e ",
            bases_mars$n_valid_interno, " VALID internos."
          )
        }
      }

      for (cenario in cenarios) {
        if (cenario == "com_MARS") {
          if (is.null(bases_mars)) {
            train_cov <- matrix(numeric(0), nrow = nrow(Xtrain), ncol = 0L,
                                dimnames = list(rownames(Xtrain), NULL))
            test_cov <- matrix(numeric(0), nrow = 1L, ncol = 0L,
                               dimnames = list(excluido, NULL))
          } else {
            train_cov <- bases_mars$train
            test_cov <- bases_mars$test
          }
        } else {
          train_cov <- Xtrain
          test_cov <- Xtest
        }
        stopifnot(identical(rownames(Lambda), rownames(train_cov)),
                  !excluido %in% rownames(train_cov))

        ncomp_max <- 0L
        pls_erro <- if (cenario == "com_MARS") erro_mars else ""
        if (ncol(train_cov) > 0L) {
          # Usar somente componentes matematicamente identificaveis no treino.
          # Se o MARS produzir <10 bases independentes, os ncomp restantes
          # ficam NA (sem acrescentar colunas originais, para preservar o metodo).
          train_pad <- scale(train_cov)
          ncomp_max <- min(nrow(train_cov) - 1L, ncol(train_cov),
                           qr(train_pad)$rank)
          ncomp_max <- min(ncomp_max, max(n.comp.pls))
        }
        componentes_disponiveis <- n.comp.pls[n.comp.pls <= ncomp_max]
        if (length(componentes_disponiveis)) {
          pls_resultado <- tryCatch({
            df_pls <- data.frame(lambda = I(Lambda), CovAmb = I(train_cov))
            mod_pls <- pls::plsr(
              lambda ~ CovAmb, data = df_pls,
              ncomp = max(componentes_disponiveis), validation = "none",
              scale = TRUE, na.action = "na.fail"
            )
            predict(mod_pls,
                    newdata = data.frame(CovAmb = I(test_cov)),
                    ncomp = componentes_disponiveis)
          }, error = function(e) e)
          if (inherits(pls_resultado, "error")) {
            pls_erro <- conditionMessage(pls_resultado)
            warning("PLS falhou FA", k, ", ", cenario, " sem ", excluido,
                    ": ", pls_erro, call. = FALSE)
            componentes_disponiveis <- integer(0)
          } else {
            pred.loads <- pls_resultado
            if (any(!is.finite(pred.loads))) {
              pls_erro <- "A PLS retornou cargas nao finitas."
              componentes_disponiveis <- integer(0)
            } else {
              for (u in seq_along(componentes_disponiveis)) {
                j <- componentes_disponiveis[u]
                pred_lambda <- matrix(pred.loads[, , u], nrow = 1L, ncol = k)
                pred_gen <- as.numeric(pred_lambda %*% t(Escores))
                gen_ids <- rownames(Escores)
                common_gen <- intersect(name.gen, gen_ids)
                if (all(is.finite(pred_gen))) {
                  pidx <- match(j, n.comp.pls)
                  pred_current[[cenario]][[pidx]][common_gen, excluido] <-
                    pred_gen[match(common_gen, gen_ids)]
                } else {
                  pls_erro <- "Predicoes geneticas nao finitas."
                  componentes_disponiveis <- setdiff(componentes_disponiveis, j)
                }
              }
            }
          }
        }
        if (length(componentes_disponiveis) < length(n.comp.pls) && verbose) {
          message("  ", cenario, " FA", k, " sem ", excluido,
                  ": componentes PLS concluidos: ",
                  if (length(componentes_disponiveis))
                    paste(componentes_disponiveis, collapse = ",") else "nenhum",
                  "; os demais permanecem NA.")
        }
        fold_diagnostics[[paste(k, excluido, cenario, sep = "|")]] <- data.frame(
          factors = k, Ambiente = excluido, cenario = cenario,
          n_train_env = nrow(train_cov), n_covariaveis = ncol(train_cov),
          n_MARS = if (cenario == "com_MARS") ncol(train_cov) else NA_integer_,
          n_complementares = if (cenario == "com_MARS") 0L else NA_integer_,
          ncomp_max = if (length(componentes_disponiveis))
            max(componentes_disponiveis) else 0L,
          max_ncomp = fold.cov[[excluido]]$limit,
          converged = isTRUE(mod$converge),
          PLS_valida = length(componentes_disponiveis) > 0L,
          MARS_aplicado = cenario == "com_MARS" && !is.null(bases_mars),
          n_train_obs = nrow(dat.train),
          n_test_obs = sum(as.character(phen.data$env) == excluido),
          test_present_in_fa = excluido %in% rownames(Lambda),
          test_present_in_pls = excluido %in% rownames(train_cov),
          erro = pls_erro,
          stringsAsFactors = FALSE
        )
      }
      checkpoint()
    }

    # Somente agora que TODAS as predicoes LOEO estao prontas: FA completo.
    reference <- NULL
    reference_model <- NULL
    erro_referencia <- ""
    if (compare.full.blups) {
      tentativa_ref <- tryCatch({
        reference_model <- fit_fa(
          phen.data, paste0("Referencia BLUP_BLUP FA", k),
          tipo_ajuste = "Referencia_completa"
        )
        ref <- as.data.frame(fa.outs.fun(
          reference_model, name.env = "env", name.gen = "gen"
        )$blups)
        if (!all(c("env", "gen", "conditional") %in% names(ref))) {
          stop("BLUP completo sem env, gen ou conditional.")
        }
        ref$env <- as.character(ref$env)
        ref$gen <- as.character(ref$gen)
        if (anyDuplicated(ref[c("env", "gen")]) ||
            !setequal(unique(ref$env), name.env)) {
          stop("BLUP completo com pares duplicados ou ambientes distintos.")
        }
        ref
      }, error = function(e) e)
      if (inherits(tentativa_ref, "error")) {
        erro_referencia <- conditionMessage(tentativa_ref)
        warning("BLUP_BLUP indisponivel FA", k, ": ", erro_referencia,
                call. = FALSE)
      } else {
        reference <- tentativa_ref
        # O FA completo foi ajustado SOMENTE apos todas as predicoes LOEO.
        # Reaproveita-lo como FA de treinamento final, sem reajuste redundante.
        pasta_fa_final <- file.path(out.dir, "modelos_finais", "FA")
        dir.create(pasta_fa_final, recursive = TRUE, showWarnings = FALSE)
        fa_extraida <- fa.outs.fun(reference_model, name.env = "env", name.gen = "gen")
        saveRDS(list(
          k = k, rot.loads = as.matrix(fa_extraida$rot.loads),
          rot.scores = as.matrix(fa_extraida$rot.scores),
          env_treino = name.env, gen_treino = name.gen,
          converged = isTRUE(reference_model$converge)
        ), file.path(pasta_fa_final, paste0("FA", k, "_extraido.rds")))
        tryCatch(
          saveRDS(reference_model,
                  file.path(pasta_fa_final, paste0("FA", k, "_ASReml.rds"))),
          error = function(e) warning(
            "Nao foi possivel serializar o ASReml completo FA", k,
            ": ", conditionMessage(e), call. = FALSE
          )
        )
      }
    }
    reference_diagnostics[[paste0("FA", k)]] <- data.frame(
      factors = k, solicitada = compare.full.blups,
      FA_convergiu = if (compare.full.blups)
        isTRUE(reference_model$converge) else NA,
      BLUPs_disponiveis = !is.null(reference),
      erro = gsub("[\r\n\t]+", " ", erro_referencia),
      stringsAsFactors = FALSE
    )
    checkpoint()

    metodos <- if (compare.full.blups) c("BLUP_BLUP", "BLUP_fen") else "BLUP_fen"
    for (cenario in cenarios) {
      for (u in seq_along(n.comp.pls)) {
        j <- n.comp.pls[u]
        predito <- pred_current[[cenario]][[u]]
        linhas_env <- list()
        for (env_id in name.env) {
          if (compare.full.blups) {
            if (is.null(reference)) {
              linhas_env[[paste0(env_id, "_blup")]] <- pair_metrics(
                env_id, "BLUP_BLUP", numeric(0), numeric(0), character(0)
              )
            } else {
              ref <- reference[reference$env == env_id, , drop = FALSE]
              p <- predito[match(ref$gen, name.gen), env_id]
              linhas_env[[paste0(env_id, "_blup")]] <- pair_metrics(
                env_id, "BLUP_BLUP", p, ref$conditional, ref$gen
              )
            }
          }
          obs <- phen.data[phen.data$env == env_id, , drop = FALSE]
          g <- as.character(obs$gen)
          p <- predito[match(g, name.gen), env_id]
          linhas_env[[paste0(env_id, "_fen")]] <- pair_metrics(
            env_id, "BLUP_fen", p, obs$prod, g
          )
        }
        por_env <- bind_df(linhas_env)
        rownames(por_env) <- NULL
        por_env$factors <- k
        por_env$ncomp <- j
        por_env$cenario <- cenario
        resumo_cv <- bind_df(lapply(metodos, function(metodo) {
          z <- por_env[por_env$method == metodo, , drop = FALSE]
          pearson <- desc_cor(z$Pearson_Correlation)
          spearman <- desc_cor(z$Spearman_Correlation)
          data.frame(
            method = metodo,
            cv.spearman = unname(spearman["mean"]),
            cv.pearson = unname(pearson["mean"]),
            sd.spearman = unname(spearman["sd"]),
            se.spearman = unname(spearman["se"]),
            sd.pearson = unname(pearson["sd"]),
            se.pearson = unname(pearson["se"]),
            n_env.spearman = unname(spearman["n"]),
            n_env.pearson = unname(pearson["n"]),
            factors = k, ncomp = j, cenario = cenario,
            reference_available = if (metodo == "BLUP_BLUP")
              !is.null(reference) else NA,
            row.names = NULL
          )
        }))
        pred_df <- data.frame(gen = name.gen, predito,
                              check.names = FALSE, row.names = NULL)
        key <- paste0(cenario, "_FA", k, "_ncomp", j)
        base <- file_base(cenario, k, j)
        write_tab(resumo_cv, paste0(base, ".txt"))
        write_tab(por_env, paste0(base, "_por_ambiente.txt"))
        write_tab(pred_df, paste0(base, "_predicoes.txt"))
        summaries[[key]] <- resumo_cv
        by_environment[[key]] <- por_env
        predictions[[key]] <- pred_df
      }
    }
    checkpoint()
  }

  tb <- convergence_tables()
  result <- list(
    summary = bind_df(summaries),
    by_environment = bind_df(by_environment),
    predictions = predictions,
    diagnostics = bind_df(fold_diagnostics),
    mars_features = bind_df(features),
    mars_resumo = bind_df(mars_resumos),
    fit_attempts = bind_df(fit_attempts),
    reference_diagnostics = bind_df(reference_diagnostics),
    convergencia = tb$detalhada,
    resumo_convergencia = tb$resumo
  )
  for (nm in c("summary", "by_environment", "diagnostics")) {
    if (!is.null(result[[nm]])) rownames(result[[nm]]) <- NULL
  }
  if (!is.null(result$mars_features)) rownames(result$mars_features) <- NULL
  if (verbose) message("LOEO concluida. Resultados: ", out.dir)
  invisible(result)
}

#==============================================================================
# 8. Executar: FA2/FA3; PLS 2:10; dois cenarios comparaveis
#==============================================================================

fatores_cv <- c(2L, 3L)
componentes_cv <- 2:10

descricao <- paste0(
  "Cafe_bienio_",
  paste(sort(unique(as.character(df_filter_cultivars$bienio))), collapse = "_"),
  "_jan15dias_LOEO_retreino_comparacao_MARS"
)
pasta_cv <- RESULT_DIR

results <- cross.val.loeo(
  factors = fatores_cv,
  phen.data = dados_cv,
  env.data = CovAmb0,
  description = descricao,
  n.comp.pls = componentes_cv,
  out.dir = pasta_cv,
  fa.outs.fun = fa.outs,
  overwrite = FALSE,
  maxit = 1000L,
  compare.full.blups = TRUE,
  cenarios = c("sem_MARS", "com_MARS"),
  max.features.por.FA = 10000L,
  max.nos.candidatos = 30L,
  min.ambientes.lado = 5L,
  min.reducao.relativa.RSS = 0.001,
  usar.interacoes = TRUE,
  prop.validacao = 0.20,
  seed.validacao = 123L
)

cv.results <- results$summary
print(cv.results)
print(results$by_environment)
print(results$diagnostics)
print(results$reference_diagnostics)
print(knitr::kable(results$resumo_convergencia))

# Auditoria LOEO para ambos os cenarios.
stopifnot(
  nrow(results$diagnostics) ==
    length(fatores_cv) * nlevels(dados_cv$env) * 2L,
  all(results$diagnostics$n_train_env == nlevels(dados_cv$env) - 1L),
  all(results$diagnostics$converged),
  # Ncomp sem bases suficientes MARS fica NA e e identificado nos diagnosticos.
  !any(results$diagnostics$test_present_in_fa),
  !any(results$diagnostics$test_present_in_pls),
  all(results$diagnostics$PLS_valida[results$diagnostics$cenario == "sem_MARS"])
)

write.table(
  results$summary, file.path(COMP_DIR, "resumo_LOEO_PLS_MARS.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)
write.table(
  results$by_environment, file.path(COMP_DIR, "LOEO_por_ambiente.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)
write.table(
  results$diagnostics, file.path(DIAG_DIR, "auditoria_treino_LOEO.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)
write.table(
  results$mars_features,
  file.path(RESULT_DIR, "com_MARS", "variaveis", "funcoes_base_MARS_por_rodada.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)
write.table(
  results$mars_resumo,
  file.path(RESULT_DIR, "com_MARS", "variaveis", "resumo_MARS_por_fator_e_rodada.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)
write.table(
  results$reference_diagnostics,
  file.path(DIAG_DIR, "diagnostico_referencia_BLUP_BLUP.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)
saveRDS(results, file.path(COMP_DIR, "results_LOEO_PLS_MARS.rds"))

# Arquivos consolidados de cada cenario e FA (alem das saidas por ncomp).
for (cenario in CENARIOS) for (k in FA_NIVEIS) {
  cenario_atual <- cenario
  fator_atual <- k
  pasta.fa <- file.path(RESULT_DIR, cenario, paste0("FA", k))
  filtro_resumo <- subset(results$summary,
                          cenario == cenario_atual & factors == fator_atual)
  filtro_env <- subset(results$by_environment,
                       cenario == cenario_atual & factors == fator_atual)
  filtro_diag <- subset(results$diagnostics,
                        cenario == cenario_atual & factors == fator_atual)
  write.table(filtro_resumo, file.path(pasta.fa, "resumo_FA.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  write.table(filtro_env, file.path(pasta.fa, "correlacoes_por_ambiente_FA.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  write.table(filtro_diag, file.path(pasta.fa, "auditoria_FA.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

# Comparacao direta com a MESMA referencia, ambientes e numero de componentes.
comparacao <- results$summary %>%
  filter(method == "BLUP_fen") %>%
  select(cenario, factors, ncomp, cv.pearson, cv.spearman,
         sd.pearson, sd.spearman, n_env.pearson) %>%
  tidyr::pivot_wider(
    names_from = cenario,
    values_from = c(cv.pearson, cv.spearman, sd.pearson, sd.spearman,
                    n_env.pearson)
  ) %>%
  mutate(
    delta_pearson = cv.pearson_com_MARS - cv.pearson_sem_MARS,
    delta_spearman = cv.spearman_com_MARS - cv.spearman_sem_MARS
  ) %>%
  arrange(factors, ncomp)
print(knitr::kable(comparacao, digits = 4))
write.table(
  comparacao, file.path(COMP_DIR, "comparacao_sem_MARS_vs_com_MARS.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)
# Grafico comparativo, sem escolher o "melhor" modelo a partir do teste LOEO.
grafico_comparativo <- results$summary %>%
  dplyr::filter(method == "BLUP_fen") %>%
  dplyr::select(cenario, factors, ncomp, cv.pearson, cv.spearman) %>%
  tidyr::pivot_longer(
    cols = c(cv.pearson, cv.spearman),
    names_to = "metrica", values_to = "correlacao"
  ) %>%
  dplyr::mutate(
    metrica = recode(metrica,
                     "cv.pearson" = "Pearson",
                     "cv.spearman" = "Spearman"),
    cenario = factor(cenario, levels = c("sem_MARS", "com_MARS"))
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = ncomp, y = correlacao,
                              color = cenario, group = cenario)) +
  ggplot2::geom_line(linewidth = 0.85) +
  ggplot2::geom_point(size = 2.1) +
  ggplot2::facet_grid(metrica ~ factors,
                      labeller = ggplot2::labeller(factors = function(x) paste0("FA", x))) +
  ggplot2::scale_x_continuous(breaks = componentes_cv) +
  ggplot2::labs(
    x = "Numero de componentes PLS", y = "Correlacao media LOEO",
    color = "Cenario", title = "GIS-FA: PLS original vs MARS.forward + PLS",
    subtitle = "FA reajustada em cada LOEO | MARS.forward com funcoes-base | BLUP_fen"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold"),
    strip.text = ggplot2::element_text(face = "bold")
  )
ggplot2::ggsave(
  file.path(COMP_DIR, "comparacao_PLS_MARS_Pearson_Spearman.png"),
  plot = grafico_comparativo, width = 12, height = 6.5, dpi = 300, bg = "white"
)
print(grafico_comparativo)

cat("Resultados salvos em:", pasta_cv, "\n")
if (any(!results$resumo_convergencia$todas_rodadas_concluidas)) {
  warning("Alguns ncomp ou rodadas MARS nao produziram predicoes. ",
          "Confira n_env.pearson, resultados por ambiente e auditoria; ",
          "nao interprete medias parciais como LOEO completa.", call. = FALSE)
}
################################################################################

cat("Cenarios separados em:", file.path(RESULT_DIR, "sem_MARS"), " e ", file.path(RESULT_DIR, "com_MARS"), "\n")

#==============================================================================
# 9. Salvar modelos finais para predicao de ambientes nao testados
#==============================================================================
# Usa FA de referencia treinada com todos ambientes, gerada SOMENTE apos LOEO.
# Nao altera as correlacoes da validacao.
script_utils <- file.path(WORK_DIR, "gisfa_funcoes_modelo.R")
if (!file.exists(script_utils)) stop("Coloque gisfa_funcoes_modelo.R no diretorio do PBS: ", WORK_DIR)
source(script_utils, local = FALSE)

selection <- escolher_configuracoes(results, nlevels(dados_cv$env))
utils::write.table(selection, file.path(RESULT_DIR, "modelos_finais", "selecao_configuracoes_LOEO.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

salvar_todos_modelos_finais(
  pasta = file.path(RESULT_DIR, "modelos_finais"),
  fa_niveis = fatores_cv,
  cenarios = CENARIOS,
  X = CovAmb0,
  selecao = selection,
  caracteristicas = caracteristicas,
  bienio = "2024-2025",
  periodos = sort(unique(as.integer(env.info_15dias$periodo))),
  max_comp_pls = max(componentes_cv),
  mars_forward = mars.forward,
  max_features_mars = 10000L,
  max_nos = 30L,
  min_lado = 5L,
  min_rss_rel = 0.001,
  prop_validacao = 0.20,
  seed_validacao = 123L
)
cat("Modelos finais exportados para: ", file.path(RESULT_DIR, "modelos_finais"), "\n")
