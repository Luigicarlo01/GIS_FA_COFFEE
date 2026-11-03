################################################################################
# GIS-FA — predicao para novos ambientes, SEM fenotipos e SEM ASReml
#
# Uso:
# Rscript --vanilla predizer_ambientes_novos.R \
#  /caminho/para/modelos_finais \
#  /caminho/para/novos_dados_diarios.RDS \
#  2026-2027 \
#  /caminho/para/resultados_predicao
#
# Usa os modelos FA2/FA3 com e sem MARS cujo ncomp foi selecionado na LOEO.
# Para inspecao de modelos nao validados, e possivel informar NCOMP_MANUAL=3
# no ambiente; ainda assim a interpretacao exige cautela.
################################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 4L) {
  stop("Uso: Rscript predizer_ambientes_novos.R ",
       "<pasta_modelos_finais> <arquivo_novos_dados.RDS> <bienio_ex_2026-2027> ",
       "[pasta_saida]", call. = FALSE)
}
modelos_dir <- path.expand(args[[1L]])
arquivo_ambiental <- path.expand(args[[2L]])
bienio_novo <- args[[3L]]
pasta_saida <- if (length(args) == 4L) path.expand(args[[4L]]) else
  file.path(dirname(modelos_dir), "predicoes_novas", paste0("bienio_", bienio_novo))
manual <- Sys.getenv("NCOMP_MANUAL", unset = "")
manual <- if (nzchar(manual)) suppressWarnings(as.integer(manual)) else NULL
if (!is.null(manual) && is.na(manual)) stop("NCOMP_MANUAL deve ser um inteiro")

for (pkg in c("dplyr", "tidyr", "pls")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Pacote ausente: ", pkg)
}
if (!dir.exists(modelos_dir)) stop("Diretorio de modelos nao encontrado: ", modelos_dir)
if (!file.exists(arquivo_ambiental)) stop("Arquivo RDS nao encontrado: ", arquivo_ambiental)
dir.create(pasta_saida, recursive = TRUE, showWarnings = FALSE)

# O diretorio com o script auxiliar faz parte do pacote de trabalho.
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
auxiliar <- file.path(script_dir, "gisfa_funcoes_modelo.R")
if (!file.exists(auxiliar)) stop("Arquivo auxiliar ausente: ", auxiliar)
source(auxiliar)

arquivos <- list.files(modelos_dir, pattern = "^modelo_final\\.rds$",
                       recursive = TRUE, full.names = TRUE)
if (!length(arquivos)) stop("Nenhum modelo_final.rds encontrado em ", modelos_dir)
modelos <- lapply(arquivos, readRDS)
nomes <- vapply(modelos, function(m) paste0(m$cenario, "_FA", m$fa), character(1L))
names(modelos) <- nomes

# Todos os modelos da execucao usam a mesma transformacao ambiental.
referencia <- modelos[[1L]]
for (m in modelos) {
  if (!identical(m$caracteristicas_diarias, referencia$caracteristicas_diarias) ||
      !identical(m$periodos, referencia$periodos)) {
    stop("Os modelos nao foram treinados com o mesmo processamento ambiental")
  }
}
cat("Lendo arquivo ambiental: ", arquivo_ambiental, "\n", sep = "")
novos_diarios <- readRDS(arquivo_ambiental)
novo_wide <- agregar_ciclos_15dias(
  novos_diarios, bienio = bienio_novo,
  caracteristicas = referencia$caracteristicas_diarias,
  periodos_esperados = referencia$periodos
)
mat <- novo_wide[, setdiff(names(novo_wide), c("env_novo", "lat", "lon")), drop = FALSE]
rownames(mat) <- novo_wide$env_novo
mat <- as.matrix(mat)
if (!is.numeric(mat)) stop("Matriz de dados ambientais agregados nao numerica")
utils::write.table(novo_wide, file.path(pasta_saida, "covariaveis_15dias_novas.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

resultados <- list()
status <- list()
for (i in seq_along(modelos)) {
  m <- modelos[[i]]
  nome <- names(modelos)[i]
  subdir <- file.path(pasta_saida, m$cenario, paste0("FA", m$fa))
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  ncomp <- if (is.null(manual)) m$ncomp_recomendado else manual
  if (is.na(ncomp) || !ncomp %in% m$ncomp_disponiveis) {
    aviso <- if (is.na(ncomp)) "Sem ncomp aprovado pela LOEO" else "NCOMP_MANUAL nao disponivel neste modelo"
    warning(nome, ": ", aviso, "; predito omitido")
    status[[nome]] <- data.frame(modelo = nome, status = "nao_predito", detalhe = aviso)
    next
  }
  if (!isTRUE(m$validacao_completa)) warning(nome,
    ": predicao com ncomp manual SEM validacao LOEO completa")
  tentativa <- tryCatch(predizer_efeitos_novos(m, mat, ncomp = ncomp),
                        error = function(e) e)
  if (inherits(tentativa, "error")) {
    aviso <- conditionMessage(tentativa)
    warning(nome, ": ", aviso)
    status[[nome]] <- data.frame(modelo = nome, status = "erro", detalhe = aviso)
    next
  }
  efeitos <- tentativa$efeitos
  pred_long <- as.data.frame(as.table(efeitos), stringsAsFactors = FALSE)
  names(pred_long) <- c("env_novo", "gen", "efeito_genetico_relativo")
  pred_long <- dplyr::left_join(pred_long,
                                novo_wide[, c("env_novo", "lat", "lon")],
                                by = "env_novo")
  pred_long <- dplyr::group_by(pred_long, env_novo)
  pred_long <- dplyr::mutate(pred_long,
    ranking = rank(-efeito_genetico_relativo, ties.method = "min"))
  pred_long <- dplyr::arrange(pred_long, env_novo, ranking, gen)
  pred_long <- dplyr::ungroup(pred_long)
  pred_long$FA <- m$fa
  pred_long$cenario <- m$cenario
  pred_long$ncomp <- ncomp
  pred_long$validacao_completa <- m$validacao_completa
  utils::write.table(pred_long, file.path(subdir, "ranking_cultivares.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  utils::write.table(
    data.frame(env_novo = rownames(tentativa$cargas), tentativa$cargas,
               row.names = NULL, check.names = FALSE),
    file.path(subdir, "cargas_fatoriais_previstas.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  saveRDS(tentativa, file.path(subdir, "predicoes_completas.rds"))
  resultados[[nome]] <- pred_long
  status[[nome]] <- data.frame(modelo = nome, status = "concluido",
                                detalhe = paste("ncomp", ncomp))
  cat("Concluido: ", nome, " (", nrow(efeitos), " ambientes, ",
      ncol(efeitos), " cultivares)\n", sep = "")
}
if (length(status)) utils::write.table(do.call(rbind, status),
  file.path(pasta_saida, "status_predicoes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
if (length(resultados)) {
  utils::write.table(do.call(rbind, resultados),
    file.path(pasta_saida, "ranking_todos_modelos.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
cat("Resultado em: ", pasta_saida, "\n", sep = "")
if (!length(resultados)) stop("Nenhum modelo gerou predicoes; veja status_predicoes.tsv")
