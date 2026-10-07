# ==============================================================================
# TÜSEB B2 Ön Çalışma - Tümörle İlişkili M2 Makrofajlar (mo-Mac / TAM)
# Hedef Gen Diferansiyel İfade Analizi ve 4 Panelli Görselleştirme
# Kohort: GSE131907 | Karşılaştırma: Tümör (tLung) vs Normal Akciğer (nLung)
# ==============================================================================

# 0. Gerekli Kütüphaneler ----------------------------------------------------
pkgs <- c("data.table", "dplyr", "tibble", "openxlsx", "ggplot2", "ggrepel", "patchwork", "scales")
new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(openxlsx)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

# 1. Dosya Yolları ve Hedef Genler -------------------------------------------
work_dir <- getwd()

anno_path      <- file.path(work_dir, "data/GSE131907_Lung_Cancer_cell_annotation.txt")
clean_rds_path <- file.path(work_dir, "data/GSE131907_matrix_unpacked.rds")
output_xlsx    <- file.path(work_dir, "M2_Makrofaj_Tumor_vs_Normal_Hedef_Genler_DE.xlsx")

target_genes <- c(
  "FABP5", "LDLR", "HMGCR", "CD36", "SQLE", "PPARGC1B", "SLC7A11", "SLC27A1",
  "FABP4", "PNPLA7", "PLIN2", "DHCR7", "HADHA", "DHCR24", "SLC27A5", "ELOVL6",
  "ACSL6", "SREBF2", "GPX4", "ACSL1", "ACAA1", "ELOVL4", "SREBF1", "HMGCS1",
  "FABP3", "HAVCR2", "SOAT1", "ACAD10", "CPT1A", "ACADS", "FITM2", "ACADM",
  "ACSL4", "LIPE", "LPAR5", "SCARB1", "FASN", "ABCA1", "ABCA3", "SLC27A4",
  "ACACA", "SCD5", "CROT", "PPARA", "PNPLA4", "PPARD", "SCD", "ABCA10", "PLIN4"
)

# 2. İfade Matrisini Yükleme -------------------------------------------------
if (!exists("expr_mat")) {
  cat("İfade matrisi yükleniyor...\n")
  expr_mat <- readRDS(clean_rds_path)
}
cat(sprintf("İfade Matrisi: %d gen x %d hücre\n", nrow(expr_mat), ncol(expr_mat)))

# 3. Anotasyon Dosyasını Okuma ve Barkod Formatını Eşleştirme ----------------
meta <- fread(anno_path, data.table = FALSE)
colnames(meta)[1] <- "Index_ID"

if (length(intersect(meta$Index_ID, colnames(expr_mat))) > 0) {
  meta$Cell_Key <- meta$Index_ID
} else if ("Barcode" %in% colnames(meta) && length(intersect(meta$Barcode, colnames(expr_mat))) > 0) {
  meta$Cell_Key <- meta$Barcode
} else {
  meta$Cell_Key <- paste0(meta$Sample, "_", meta$Barcode)
}

# 4. M2 / TAM Hücrelerini Seçme ----------------------------------------------
# GSE131907'de M2 profili sergileyen primer kitle mo-Mac (monosit kökenli makrofaj) popülasyonudur.
m2_types <- c("mo-Mac") 

sub_meta <- meta %>%
  filter(
    Cell_subtype %in% m2_types,
    Sample_Origin %in% c("tLung", "nLung"),
    Cell_Key %in% colnames(expr_mat)
  ) %>%
  mutate(Group = ifelse(Sample_Origin == "tLung", "Tümör", "Normal"))

n_tumor <- sum(sub_meta$Group == "Tümör")
n_normal <- sum(sub_meta$Group == "Normal")
cat(sprintf("\nAnalize giren Hücre Sayıları -> Tümör: %d | Normal: %d\n", n_tumor, n_normal))

# 5. Diferansiyel Ekspresyon & Oran Hesaplamaları ----------------------------
common_genes <- intersect(target_genes, rownames(expr_mat))
missing_genes <- setdiff(target_genes, rownames(expr_mat))

if (length(missing_genes) > 0) {
  cat("Matriste yer almayan genler:\n", paste(missing_genes, collapse = ", "), "\n")
}

sub_expr <- expr_mat[common_genes, sub_meta$Cell_Key, drop = FALSE]

tumor_cells <- sub_meta$Cell_Key[sub_meta$Group == "Tümör"]
normal_cells <- sub_meta$Cell_Key[sub_meta$Group == "Normal"]

cat("\nWilcoxon testleri ve ekspresyon yüzdeleri hesaplanıyor...\n")

de_results <- lapply(common_genes, function(gene) {
  t_vals <- as.numeric(sub_expr[gene, tumor_cells])
  n_vals <- as.numeric(sub_expr[gene, normal_cells])
  
  mean_t <- mean(t_vals)
  mean_n <- mean(n_vals)
  
  pct_1 <- mean(t_vals > 0)
  pct_2 <- mean(n_vals > 0)
  
  w_test <- tryCatch(
    wilcox.test(t_vals, n_vals, exact = FALSE),
    error = function(e) list(p.value = NA)
  )
  
  data.frame(
    Gen = gene,
    Ortalama_Tumor = round(mean_t, 4),
    Ortalama_Normal = round(mean_n, 4),
    Log2FC = round(mean_t - mean_n, 4),
    Oran_Tumor = round(pct_1, 4),
    Oran_Normal = round(pct_2, 4),
    Oran_Fark = round(pct_1 - pct_2, 4),
    P_Degeri = w_test$p.value,
    stringsAsFactors = FALSE
  )
}) %>% 
  bind_rows() %>%
  mutate(
    Duzeltilmis_P = p.adjust(P_Degeri, method = "BH"),
    Duzenlenme = case_when(
      Log2FC > 0.20 & Duzeltilmis_P < 0.05 ~ "Tümörde Artan",
      Log2FC < -0.20 & Duzeltilmis_P < 0.05 ~ "Tümörde Azalan",
      TRUE ~ "Anlamsız"
    )
  ) %>%
  arrange(Duzeltilmis_P, desc(abs(Log2FC)))

# 6. Excel Olarak Kaydetme ----------------------------------------------------
cat("\nExcel dosyası yazılıyor...\n")
wb <- createWorkbook()
sheet_name <- "M2_Makrofaj_DE_Sonuclari"
addWorksheet(wb, sheet_name)

header_style <- createStyle(
  fontSize = 11, fontColour = "#FFFFFF", fgFill = "#1F4E78",
  textDecoration = "bold", halign = "center", valign = "center",
  border = "TopBottomLeftRight"
)
num_style <- createStyle(numFmt = "0.0000", halign = "right")
sci_style <- createStyle(numFmt = "0.00E+00", halign = "right")

writeData(wb, sheet_name, de_results, startRow = 1, startCol = 1, headerStyle = header_style)
addStyle(wb, sheet_name, style = num_style, rows = 2:(nrow(de_results) + 1), cols = 2:7, gridExpand = TRUE)
addStyle(wb, sheet_name, style = sci_style, rows = 2:(nrow(de_results) + 1), cols = 8:9, gridExpand = TRUE)
setColWidths(wb, sheet_name, cols = 1:ncol(de_results), widths = "auto")

saveWorkbook(wb, output_xlsx, overwrite = TRUE)
cat(sprintf("Excel tablosu kaydedildi:\n%s\n", output_xlsx))

# 7. Türkçe 4'lü Panel Figür Oluşturma ----------------------------------------
theme_nature <- function() {
  theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = element_text(size = 9, color = "grey30"),
      axis.title = element_text(face = "bold", size = 10),
      axis.text = element_text(color = "black", size = 9),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
    )
}

df_plot <- de_results %>%
  mutate(
    P_val_clean = ifelse(P_Degeri == 0 | is.na(P_Degeri), 1e-300, P_Degeri),
    Adj_P_val_clean = ifelse(Duzeltilmis_P == 0 | is.na(Duzeltilmis_P), 1e-300, Duzeltilmis_P),
    neg_log10_padj = -log10(Adj_P_val_clean),
    Anlamlilik = factor(case_when(
      Log2FC >= 0.20 & Duzeltilmis_P < 0.05 ~ "Tümörde Artan",
      Log2FC <= -0.10 & Duzeltilmis_P < 0.05 ~ "Tümörde Azalan",
      TRUE ~ "Anlamsız"
    ), levels = c("Tümörde Artan", "Tümörde Azalan", "Anlamsız"))
  )

label_genes <- df_plot %>%
  filter(abs(Log2FC) > 0.15 | Duzeltilmis_P < 1e-5 | Gen %in% c("CD36", "FABP4", "FABP5", "PLIN2", "PPARG", "ABCA1", "CPT1A", "GPX4", "FASN")) %>%
  pull(Gen) %>%
  unique()
label_genes <- head(label_genes, 12)

# PANEL A: Yanardağ (Volcano) Grafiği
pA <- ggplot(df_plot, aes(x = Log2FC, y = neg_log10_padj)) +
  geom_point(aes(color = Anlamlilik, size = Oran_Tumor), alpha = 0.8) +
  scale_color_manual(values = c(
    "Tümörde Artan" = "#B22222",
    "Tümörde Azalan" = "#1F78B4",
    "Anlamsız" = "grey70"
  )) +
  scale_size_continuous(range = c(1.5, 4.5), labels = scales::percent_format(accuracy = 1)) +
  geom_vline(xintercept = c(-0.20, 0.20), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_text_repel(
    data = filter(df_plot, Gen %in% label_genes),
    aes(label = Gen),
    box.padding = 0.4,
    max.overlaps = 20,
    fontface = "bold.italic",
    size = 3.3
  ) +
  labs(
    title = "A. M2 Makrofaj (TAM) Hedef Gen Yanardağ Dağılımı",
    subtitle = "Volcano Dağılımı: Tümör ve Normal Akciğer Karşılaştırması",
    x = "Katlanma Değişimi (Log2FC: Tümör / Normal)",
    y = "-Log10 Düzeltilmiş P-Değeri",
    color = "Ekspresyon Durumu",
    size = "Tümörde İfade %"
  ) +
  theme_nature()

# PANEL B: Log2FC Çubuk Grafiği
top_b <- df_plot %>%
  filter(abs(Log2FC) > 0.05 | Gen %in% label_genes) %>%
  arrange(desc(Log2FC)) %>%
  slice(c(1:min(10, n()), (max(1, n() - 3)):n())) %>%
  distinct(Gen, .keep_all = TRUE)

pB <- ggplot(top_b, aes(x = reorder(Gen, Log2FC), y = Log2FC, fill = Log2FC > 0)) +
  geom_col(width = 0.72, show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#D95F02", "FALSE" = "#7570B3")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_text(
    aes(
      label = sprintf("%.2f", Log2FC),
      hjust = ifelse(Log2FC > 0, -0.15, 1.15)
    ),
    size = 2.8,
    fontface = "bold"
  ) +
  labs(
    title = "B. Gen İfade Değişim Büyüklüğü",
    subtitle = "En çok değişen lipit ve metabolik aday genler",
    x = "",
    y = "Log2 Katlanma Değişimi (Log2FC)"
  ) +
  expand_limits(y = c(min(top_b$Log2FC) * 1.35, max(top_b$Log2FC) * 1.25)) +
  theme_nature()

# PANEL C: Nokta (Dot Plot) Grafiği
dot_genes <- head(label_genes, 10)
t_part <- df_plot %>% filter(Gen %in% dot_genes) %>% transmute(Gen, Grup = "Tümör", Ifade = Ortalama_Tumor, Yuzde = Oran_Tumor * 100)
n_part <- df_plot %>% filter(Gen %in% dot_genes) %>% transmute(Gen, Grup = "Normal", Ifade = Ortalama_Normal, Yuzde = Oran_Normal * 100)
dot_df <- bind_rows(t_part, n_part) %>%
  mutate(
    Gen = factor(Gen, levels = rev(dot_genes)),
    Grup = factor(Grup, levels = c("Normal", "Tümör"))
  )

pC <- ggplot(dot_df, aes(x = Grup, y = Gen)) +
  geom_point(aes(size = Yuzde, color = Ifade)) +
  scale_size_continuous(range = c(1.5, 6), breaks = c(10, 30, 60, 90), name = "Pozitif Hücre (%)") +
  scale_color_gradient(low = "#4575B4", high = "#D73027", name = "Ortalama İfade") +
  labs(
    title = "C. Hücresel İfade ve Yaygınlık",
    subtitle = "Kilit metabolik genlerde Tümör ve Normal doku kıyaslaması",
    x = "",
    y = "Hedef Genler"
  ) +
  theme_nature()

# PANEL D: Delta-Prevalans vs Log2FC Grafiği
pD <- ggplot(df_plot, aes(x = Oran_Fark * 100, y = Log2FC)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  geom_point(aes(color = Anlamlilik, size = Ortalama_Tumor), alpha = 0.85) +
  scale_color_manual(values = c(
    "Tümörde Artan" = "#B22222",
    "Tümörde Azalan" = "#1F78B4",
    "Anlamsız" = "grey70"
  )) +
  scale_size_continuous(range = c(2, 5), name = "Tümör İfade Seviyesi") +
  geom_text_repel(
    data = filter(df_plot, Gen %in% label_genes),
    aes(label = Gen),
    box.padding = 0.35,
    max.overlaps = 15,
    fontface = "bold.italic",
    size = 3.2
  ) +
  labs(
    title = "D. Fenotipik Yaygınlık ve İfade İlişkisi",
    subtitle = "Hücre frekans değişimi (Δ%) ile ekspresyon seviyesi korelasyonu",
    x = "Hücre Oranı Farkı: Tümör - Normal (Δ%)",
    y = "Log2 Katlanma Değişimi (Log2FC)"
  ) +
  theme_nature()

# Figür Birleştirme ve Kaydetme
fig <- (pA + pB) / (pC + pD) +
  plot_annotation(
    title = "KHDAK Tümörle İlişkili Makrofajlarda (mo-Mac / TAM) Lipit ve Metabolik Yeniden Programlanma",
    subtitle = "Kohort: GSE131907 | Karşılaştırma: Primer Tümör (tLung) ve Eşlenik Normal Doku (nLung)",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey35")
    )
  )

png_out <- file.path(work_dir, "M2_TAM_4Panel_DE_Figuru.png")
pdf_out <- file.path(work_dir, "M2_TAM_4Panel_DE_Figuru.pdf")

ggsave(png_out, plot = fig, width = 14, height = 12, dpi = 300)
ggsave(pdf_out, plot = fig, width = 14, height = 12, device = "pdf")
cat(sprintf("M2/TAM figürleri kaydedildi:\n1) %s\n2) %s\n", png_out, pdf_out))