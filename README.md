# TÜSEB B2 Projesi Ön Çalışma Raporu: Tek Hücre RNA-seq Analizi (GSE131907)

Bu depo, TÜSEB B2 proje başvurusu kapsamında yürütülen ön çalışma niteliğindeki biyoinformatik analiz betiklerini içermektedir. Küçük Hücreli Dışı Akciğer Kanseri (KHDAK / NSCLC) mikroyapısında tükenmiş CD8 T lenfositleri (CD8Tex) ve tümörle ilişkili makrofajlar (mo-Mac / TAM) üzerinde hedef lipit ve metabolik yolak regülatör genlerinin diferansiyel ifade analizleri gerçekleştirilmiştir.

---

## 📁 Depo İçeriği ve Analiz Boru Hattı

1. **`01_CD8Tex_Hedef_Gen_DE_Analizi.R`**
   - CD8Tex popülasyonunda tümör (`tLung`) ve tümör komşusu normal (`nLung`) dokuların karşılaştırılması.
   - Wilcoxon Rank Sum testi ve Benjamini-Hochberg (FDR) düzeltmesi.
   - 4 Panelli görselleştirme: Volcano grafiği, Log2FC büyüklüğü, Dot plot (prevalans/ifade) ve Fenotipik penetrans grafiği.

2. **`02_M2_Makrofaj_Hedef_Gen_DE_Analizi.R`**
   - Monosit kökenli M2 benzeri makrofaj (`mo-Mac`) alt kümesinde hedef metabolik genlerin taranması.
   - Hücresel prevalans farkları (Delta) ve katlanma değişimi (Log_2FC) analizleri.
   - Yüksek çözünürlüklü 4 panelli karşılaştırma figürü üretimi.

---

## 🛠️ Gereksinimler (R Paketleri)
Analizlerin tekrarlanabilirliği için aşağıdaki R kütüphaneleri kullanılmaktadır:
- `data.table`, `dplyr`, `tibble`
- `openxlsx` (Tablo çıktıları için)
- `ggplot2`, `ggrepel`, `patchwork`, `scales` (Görselleştirme için)

---

## 📊 Veri Kaynağı
- **Veri Seti:** GEO Veri Tabanı - [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907)
- **Tür:** Tek hücre RNA dizileme (scRNA-seq), 10x Genomics
# TÜSEB B2 Projesi Ön Çalışma: Tek Hücre RNA-seq Analiz Kodları

> 🔗 **Kaynak Kodlar:** [GitHub Deposundaki R Kodlarını Görüntülemek İçin Tıklayın](https://github.com/durualtinel0111/On-Calisma-scRNAseq)

### Analiz Betikleri:
- 📄 [01_CD8Tex_Hedef_Gen_DE_Analizi.R](https://github.com/durualtinel0111/On-Calisma-scRNAseq/blob/main/01_CD8Tex_Hedef_Gen_DE_Analizi.R)
- 📄 [02_M2_Makrofaj_Hedef_Gen_DE_Analizi.R](https://github.com/durualtinel0111/On-Calisma-scRNAseq/blob/main/02_M2_Makrofaj_Hedef_Gen_DE_Analizi.R)

---
