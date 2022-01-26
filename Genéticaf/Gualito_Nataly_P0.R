library("msa")
library(BiocGenerics)
library(BiocManager)
library(Biostrings)
#cargue las librerias

Estrella<-readRNAStringSet("C:/Users/nat_g/OneDrive/Escritorio/Genéticaf/first (1).fasta")
Estrella
#Ahora cargue el archivo FASTA


#PRIMER EJERCICIO
#Counting DNA Nucleotides

C2<-letterFrequency(Estrella,"C")#En realidad, este sirve para contar las
#C que hay en cada secuencia y pensaba hacerlo asi
C2

#Pero el ejercicio de Rosalind me dice que requiere de los resultados
#de conteo de las letras de esta manera
#En este caso es de la secuencia0
letterFrequency(Estrella[[1]], letters="ACGU", OR=0)

#Segunda secuencia o secuencia1
letterFrequency(Estrella[[2]], letters="ACGU", OR=0)

#Tercer secuencia o secuencia2
letterFrequency(Estrella[[3]], letters="ACGU", OR=0)

#Cuarta secuencia o secuencia3
letterFrequency(Estrella[[4]], letters="ACGU", OR=0)

#Quinta secuencia o secuencia4
letterFrequency(Estrella[[5]], letters="ACGU", OR=0)



#SEGUNDO EJERCICIO
#Computing GC Content

cont<-letterFrequency(Estrella, "GC", as.prob = TRUE)
cont #En este caso, no estoy contando la cantidad de G o C
#si no el contenido de GC que hay en cada secuencia




############
#Tengo entendido que debemos hacerlo sin librerías, porlo que asimilo que es 
#como de forma manual o como en ciclos, pero no estoy segura

estrella<-read.delim("C:/Users/nat_g/OneDrive/Escritorio/Genéticaf/tareasec.txt")
estrella 
#Decidi convertirlo a texto, para que me sea mas facil manejarlo

#Selecciono la fila 1, segun el documento que estoy utilizando

#Secuencia0
sec0<-estrella[1,]
sec0

#Secuencia1
sec1<-estrella[3,]
sec1

#Secuencia2
sec2<-estrella[5,]
sec2

#Secuencia3
sec3<-estrella[7,]
sec3

#Secuencia4
sec4<-estrella[9,]
sec4

#PRIMER EJERCICIO
#Counting DNA Nucleotides

#me cuenta cuantas "A", "C", "G" y "U" tengo en la primer secuencia
lengths(regmatches(sec0, gregexpr("A",sec0)))
lengths(regmatches(sec0, gregexpr("C",sec0)))
lengths(regmatches(sec0, gregexpr("G",sec0)))
lengths(regmatches(sec0, gregexpr("U",sec0)))

#segunda secuencia
lengths(regmatches(sec1, gregexpr("A",sec1)))
lengths(regmatches(sec1, gregexpr("C",sec1)))
lengths(regmatches(sec1, gregexpr("G",sec1)))
lengths(regmatches(sec1, gregexpr("U",sec1)))

#tercer secuencia
lengths(regmatches(sec2, gregexpr("A",sec2)))
lengths(regmatches(sec2, gregexpr("C",sec2)))
lengths(regmatches(sec2, gregexpr("G",sec2)))
lengths(regmatches(sec2, gregexpr("U",sec2)))

#cuarta secuencia
lengths(regmatches(sec3, gregexpr("A",sec3)))
lengths(regmatches(sec3, gregexpr("C",sec3)))
lengths(regmatches(sec3, gregexpr("G",sec3)))
lengths(regmatches(sec3, gregexpr("U",sec3)))

#quinta secuencia
lengths(regmatches(sec4, gregexpr("A",sec4)))
lengths(regmatches(sec4, gregexpr("C",sec4)))
lengths(regmatches(sec4, gregexpr("G",sec4)))
lengths(regmatches(sec4, gregexpr("U",sec4)))


#SEGUNDO EJERCICIO
#Computing GC Content


sec0.1<-c('U','A','A','U','A','A','U','A','A','U','A','A','U','A','A','U','A','A','U','A','A')
sec0.1 #Primero uso la funcion para concatenar mi objeto

#calculo la frecuencia relativa de la secuencia0
table(sec0.1)/length(sec0.1)


#secuencia1
sec1.1<-c('C','A','U','G','C','U','C','C','U','C','C','C','U','A','U')
sec1.1

#Frecuencia relativa
table(sec1.1)/length(sec1.1)


#secuencia2
sec2.1<-c('A','A','C','G','A','G','U','G','G')
sec2.1

#frecuencia relativa
table(sec2.1)/length(sec2.1)


#secuencia3
sec3.1<-c('U','A','C','G','A','G','G','C','G','A','G','G')

#frecuencia relativa
table(sec3.1)/length(sec3.1)


#secuencia4
sec4.1<-c('U','A','A','U','A','A','U','A','A','U','A','A','U','A','A','U','A','A','U','A','A')

#frecuencia relativa
table(sec4.1)/length(sec4.1)
