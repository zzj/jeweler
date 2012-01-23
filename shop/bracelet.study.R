source('pileup.plot.R')

a='CUFF.2721'
b='CUFF.19465'
a='CUFF.2477'
b='CUFF.17814'

pdf('result.plot.pdf')
a='CUFF.9'
b='CUFF.20306'
a='CUFF.484'
b='CUFF.9622'

gene.pileup.plot(paste('../result/merged_list/jeweler//HF_0128_M_merged/',a,'/',sep=""),
                 paste('../result/merged_list/jeweler/HF_0128_M_merged/',a,'/',a,'.landscape.plot.meta',sep=""),
                 paste('../result/merged_list/jeweler/HF_0128_M_merged/',a,'/',a,'.mismacher',sep=""),
                 a)
gene.pileup.plot(paste('../result/merged_list/jeweler//HF_0128_M_merged/',b,'/',sep=""),
                 paste('../result/merged_list/jeweler/HF_0128_M_merged/',b,'/',b,'.landscape.plot.meta',sep=""),
                 paste('../result/merged_list/jeweler/HF_0128_M_merged/',b,'/',b,'.mismacher',sep=""),
                 b)
dev.off()

a='CUFF.2414'
b='CUFF.24085'
a='CUFF.1091'
b='CUFF.24090'

a='CUFF.1091'
postscript(paste(a,'.eps',sep=""))
gene.pileup.plot(paste('../result/merged_list/jeweler//HF_0128_M_merged/',a,'/',sep=""),
                 paste('../result/merged_list/jeweler/HF_0128_M_merged/',a,'/',a,'.landscape.plot.meta',sep=""),
                 paste('../result/merged_list/jeweler/HF_0128_M_merged/',a,'/',a,'.mismacher',sep=""),
                 a)
dev.off()
b='CUFF.24090'
postscript(paste(b,'.eps',sep=""))
gene.pileup.plot(paste('../result/merged_list/jeweler//HF_0128_M_merged/',b,'/',sep=""),
                 paste('../result/merged_list/jeweler/HF_0128_M_merged/',b,'/',b,'.landscape.plot.meta',sep=""),
                 paste('../result/merged_list/jeweler/HF_0128_M_merged/',b,'/',b,'.mismacher',sep=""),
                 b)

dev.off()
