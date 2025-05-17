# auteur: Tidjani D.
# date: 5 mars 2024


#################### Fichier de départ du programme ##########################

adn = "TCGACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGACG\
GCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGAGTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGAATGCCAGCCAGC\
CAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGAACACTCTTCGACAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGA\
ACTCGACACTCTTCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCATCCCAGCGATACCC\
AGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAG\
CCAGCGAACTCGTCTGCGTTCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGC\
GATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGTATGCCAGCC\
AGCATCCCAGCGA"

codons_aa = {
    "UUU": "Phénylalanine",
    "UUC": "Phénylalanine",
    "UUA": "Leucine",
    "UUG": "Leucine",
    "CUU": "Leucine",
    "CUC": "Leucine",
    "CUA": "Leucine",
    "CUG": "Leucine",
    "AUU": "Isoleucine",
    "AUC": "Isoleucine",
    "AUA": "Isoleucine",
    "AUG": "Méthionine (Start)",
    "GUU": "Valine",
    "GUC": "Valine",
    "GUA": "Valine",
    "GUG": "Valine",
    "UCU": "Sérine",
    "UCC": "Sérine",
    "UCA": "Sérine",
    "UCG": "Sérine",
    "CCU": "Proline",
    "CCC": "Proline",
    "CCA": "Proline",
    "CCG": "Proline",
    "ACU": "Thrénine",
    "ACC": "Thrénine",
    "ACA": "Thrénine",
    "ACG": "Thrénine",
    "GCU": "Alanine",
    "GCC": "Alanine",
    "GCA": "Alanine",
    "GCG": "Alanine",
    "UAU": "Tyrosine",
    "UAC": "Tyrosine",
    "UAA": "Stop",
    "UAG": "Stop",
    "CAU": "Histidine",
    "CAC": "Histidine",
    "CAA": "Glutamine",
    "CAG": "Glutamine",
    "AAU": "Asparagine",
    "AAC": "Asparagine",
    "AAA": "Lysine",
    "AAG": "Lysine",
    "GAU": "Aspartate",
    "GAC": "Aspartate",
    "GAA": "Glutamate",
    "GAG": "Glutamate",
    "UGU": "Cystéine",
    "UGC": "Cystéine",
    "UGA": "Stop",
    "UGG": "Tryptophane",
    "CGU": "Arginine",
    "CGC": "Arginine",
    "CGA": "Arginine",
    "CGG": "Arginine",
    "AGU": "Sérine",
    "AGC": "Sérine",
    "AGA": "Arginine",
    "AGG": "Arginine",
    "GGU": "Glycine",
    "GGC": "Glycine",
    "GGA": "Glycine",
    "GGG": "Glycine"}

lettreAa = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "UGU": "C",
    "UGC": "C",
    "UGA": "*",
    "UGG": "W",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G"
}

############### Développement des fonctions et procédures #####################

### Développement de la fonction antisens(brinAdn)

# tâche : Part de brin ADN fourni et renvoie le brin d’ADN complémentaire.


def antisens(brinAdn):
    brinAdnComplementaire = []
    for i in brinAdn[::-1]: # on inverse le sens du brin pour obtenir le brin
        if i == 'A':        # complémentaire, tout en remplaçant les nucléoti-
            i = 'T'         # des par leurs équivalents.
            brinAdnComplementaire.append(i)
        elif i =='C':
            i = 'G'
            brinAdnComplementaire.append(i)
        elif i == 'T':
            i = 'A'
            brinAdnComplementaire.append(i)
        elif i == 'G':
            i = 'C'
            brinAdnComplementaire.append(i)
    return (''.join(brinAdnComplementaire))

def testAntiSens():
    assert antisens('ATCAATTGCAAGTTCC') == 'GGAACTTGCAATTGAT'
    assert antisens('GCGATCGATCGATCGA') == 'TCGATCGATCGATCGC'
    assert antisens('TCGATCGATCGATCGA') == 'TCGATCGATCGATCGA'
    assert antisens('AAGCTTGGCCAAAGCT') == 'AGCTTTGGCCAAGCTT'
    assert antisens('GATCGATCGATCGATC') == 'GATCGATCGATCGATC'
    

    
### Développement de la fonction trouveDebut(brinAdn)

# tâche : Recherche tous les codons de départ sur un brin d’ADN et renvoie un 
# tableau contenant les positions du premier nucléotide de chacun des codons.

def trouveDebut(brinAdn):          # un gène débute tout le temps par la sé-
    positionsDebut = []            # quence 'TAC'. Donc on parcourt le brin 
    for j in range(len(brinAdn)):  # d'ADN à la recherche des positions des 'T'
        if brinAdn[j:j+3] == 'TAC':# suivis de 'AC'
            positionsDebut.append(j)
    return positionsDebut


def testTrouveDebut():
    assert trouveDebut('ACTGATACGCTACGATACGATACGCTACGATACGCTACGATACGCTACGATACG\
    CTACG') == [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 59]
    assert trouveDebut('GATACGCTAGCTACGATACGCTAGCTACGATACGCTAGCTACGATACGCTAGCT\
    ACGATA') == [2, 11, 16, 25, 30, 39, 44]
    assert trouveDebut('CGTACTAGCTACGATACGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG\
    CTAGC') == [2, 9, 14]
    assert trouveDebut('ACGTACGTGACTAGCTAGCTAGCATCGATCGATCGATCGTAGCTAGCTAGCTAG\
    CTAGC') == [3]
    assert trouveDebut('CGATCGTAGCGATCGATCGTAGCGATCGTAGCGATCGTAGCGATCGTAGCGATC\
    GTAGC') == []
    


### Développement de la fonction trouveFin(brinAdn)

# tâche : Même chose que la fonction précédente mais renvoie un tableau avec 
# les positions de tous les codons de terminaison

def trouveFin(brinAdn):                     # même chose que précedemment sauf
    positionsFin = []                       # que là, on parcours le brin 
    terminaisons = ['ATT','ATC','ACT']      # d'ADN à la recherche des posi-
    for k in range(len(brinAdn)):           # tions des 'A' suivis de 'TT',
        if brinAdn[k:k + 3] in terminaisons:# 'TC' ou alors 'CT'.
            positionsFin.append(k+3)
    return positionsFin

def testTrouveFin():
    assert trouveFin('TCGACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCA\
    GCGAAGCCAGCCAGCCGA') == [6, 12]
    assert trouveFin('GCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGAGTGCCAGCCAGCCAGCCAGCGA\
    ACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGAATGCCAGCCAGC') == [63, 69]
    assert trouveFin('CAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGAACACTCTTCG\
    ACAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGA') == [28, 51, 89]
    assert trouveFin('AGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACT\
    GCGATCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAG') == [33, 56, 66, 94]
    assert trouveFin('CGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAGC\
    GTAG') == []



### Développement de la fonction trouveGene(debut, fin)

# tâche : Prend en paramètre un tableau contenant les positions de tous les 
# codons de départ et un autre tableau contenant les positions de tous les 
# codons de terminaison pour un brin d’ADN et renvoie un tableau de tuples 
# contenant la liste des gènes (début et fin) trouvés sur un brin.
 
def trouveGene(positionsDebut,positionsFin): # on vérifie grâce à un modulo
    positionsGenes = []                      # si le nombre de nucléotides
    for m in range(len(positionsDebut)):     # entre une position de début 
        for n in range(len(positionsFin)):   # et celle d'une fin est un multi-
                                             #ple de 3.
            if positionsFin[n] > positionsDebut[m] and \
            (positionsFin[n] - positionsDebut[m]) % 3 == 0:
                positionsGenes.append((positionsDebut[m], positionsFin[n]))
                break
    return positionsGenes

def testTrouveGene():
    assert trouveGene([38, 74, 402],[3, 9, 154, 160, 226, 249, 283, 304, 311, 
                                     343, 379, 392, 437, 460, 466, 494, 517, 
                                     556, 592, 720]) == [(38, 311), (74, 311),
                                                         (402, 720)]
    assert trouveGene([9, 14],[2, 54]) == [(9, 54)]
    assert trouveGene([6, 18, 31, 45, 57],
                      [12, 25, 38, 51, 63]) == [(6, 12), (18, 51), 
                                                (45, 51), (57, 63)]
    assert trouveGene([0, 12, 24],[9, 21, 33]) == [(0, 9), (12, 21), (24, 33)]
    assert trouveGene([3, 12, 21, 30],[9, 18, 27, 36]) == [(3, 9), (12, 18), 
                                                           (21, 27), (30, 36)]




### Fonction transcrire(brinAdn)

## Développement de la fontion intermédiaire lesGenes(positionsGenes,brinAdn)

# tâche : Prend en paramètre le tableau de tuples contenant la liste des gènes
# (début et fin) trouvés sur un brin et le brin lui-même puis renvoie un 
# tableau contenant les sous-chaines de caractère du brin d’ADN.

def lesGenes(positionsGenes,brinAdn):
    mesGenes = []
    mesIndex = []
    for i in positionsGenes:
        mesIndex.append(i[0]) # mesIndex représente un tableau contenant 
        mesIndex.append(i[1]) # les positions de début et de fin des gènes
    monBooleen = False        # successivement, afin de faciliter l'extraction
    for j in mesIndex:        # des sous-chaînes de caractères correspondant
        if monBooleen:        # aux gènes.
            fin = j+1
            mesGenes.append(brinAdn[debut:fin])
            monBooleen = False
        else:                # à l'aide d'un booléen, on récupère nos gènes.
            debut = j
            monBooleen = True
    return mesGenes



## Développement de la fonction principale transcrire(brinAdn)

# tâche : Prend en paramètre la sous-chaine de caractère du brin d’ADN débutant
# au début du gène et se terminant à la fin du gène et renvoie le brin d’ARN 
# correspondant sous forme d’une chaine de caractères.

def transcrire(brinAdn):     # là où il devait y avoir des 'T', on les remplace
    monARN = []              # par des 'U'.
    for i in brinAdn:
        if i == 'A':
            i = 'U'
            monARN.append(i)
        elif i == 'C':
            i = 'G'
            monARN.append(i)
        elif i == 'T':
            i = 'A'
            monARN.append(i)
        elif i == 'G':
            i = 'C'
            monARN.append(i)
    monARN = ''.join(monARN)
    return monARN

def testTranscrire():
    assert transcrire('ACTGATACGCTACGATACGATACGCTACGAT') == \
    'UGACUAUGCGAUGCUAUGCUAUGCGAUGCUA'
    assert transcrire('GATACGCTAGCTACGATACGCTAGCTACGAT') == \
    'CUAUGCGAUCGAUGCUAUGCGAUCGAUGCUA'
    assert transcrire('AGCTACGATACGCTAGCTAGCTA') == 'UCGAUGCUAUGCGAUCGAUCGAU'
    assert transcrire('TCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGC') == \
    'AGCUAGCUAGCUAGCAUCGAUCGAUCGAUCGAUCG'
    assert transcrire('GATCGTAGCGATCGTAGCGATCGT') == 'CUAGCAUCGCUAGCAUCGCUAGCA'


### Développement de la procédure traduire(brinArn)

# tâche : Prend en paramètre un brin d’ARN (chaine de caractères) et affiche la
# protéine sous forme d’une chaine de caractères et la dessine à l’aide 
# de la tortue.

## Développement de la fonction fromArnToProteine(brinArn)

# tâche : Prend en paramètre un brin d’ARN (chaine de caractères) et affiche 
# la protéine sous forme d’une chaine de caractères.

def fromArnToProteine(brinArn):
    proteine = []
    for i in range(0,len(brinArn),3):
        if brinArn[i:i + 3] in codons_aa:
            proteine.append(codons_aa.get(brinArn[i:i + 3]))
    lastElement = proteine.pop() # pour éliminer le codon STOP.
    proteine = '-'.join(proteine)
    return proteine

def testfromArnToProteine():
    assert fromArnToProteine('AUGUAUGAUUAUUAA') == \
    'Méthionine (Start)-Tyrosine-Aspartate-Tyrosine'
    assert fromArnToProteine('AUGGGUGGGUCAUAG') == \
    'Méthionine (Start)-Glycine-Glycine-Sérine'
    assert fromArnToProteine('AUGGGUUUUGAAUAG') == \
    'Méthionine (Start)-Glycine-Phénylalanine-Glutamate'


## Développement de la procédure de dessin de la protéine 
# lettresDesCodons(brinArn)

# tâche : Prend en paramètre la protéine sous forme d’une chaine de caractères
# puis la dessine à l'aide de la Tortue.


def lettresDesCodons(brinArn):
    proteine = []
    codonsStop = ['UAA', 'UAG', 'UGA'] 
    sum = 0
    for i in range(0,len(brinArn)-1,3):
        if brinArn[i:i + 3] not in codonsStop and sum < 15:
            lettre = lettreAa.get(brinArn[i:i + 3])
            write(lettre)
            pu()
            fd(20)
            pd()
            sum += 1
        elif sum%15 == 0: # dans ce cas, on met la tortue à la ligne et on
            pu()          # réinitialise sum afin de commencer une nouvelle
            rt(90)        # ligne.
            fd(20)
            lt(90)
            bk(300)
            pd()
            sum = 0

            lettre = lettreAa.get(brinArn[i:i + 3])
            write(lettre)
            pu()
            fd(20)
            pd()
            sum += 1
            
    rt(90)
       


## Développement de la procédure de dessin carre(longueur, nombre)

# tâche : Prend deux entiers en paramètre (taille du côté du carré et l’indice 
# du carré à dessiner) et trace un carré à l’aide de la tortue.

def carre(longueur,nombre):
    temp1 = nombre // 15   # pour définir le nombre de ligne.
    temp2 = nombre % 15    # pour définir le nombre de carrés restants 
    for j in range(temp1): # à dessiner après avoir terminer les lignes.
        for i in range(15):
            fd(longueur)
            rt(90)
            fd(longueur)
            rt(90)
            fd(longueur)
            rt(90)
            fd(longueur)
            pu()
            rt(90)
            fd(longueur)
            pd()
        if not nombre<15 :
            pu()
            rt(90)
            fd(longueur)
            rt(90)
            fd(15*longueur)
            rt(180)
            pd()
    if temp2 != 0:             # pour éviter un bug lorsqu'il y'a pas de carrés 
        for k in range(temp2): # restants à dessiner après avoir compléter 
            fd(longueur)       # les lignes de carrés.
            rt(90)
            fd(longueur)
            rt(90)
            fd(longueur)
            rt(90)
            fd(longueur)
            pu()
            rt(90)
            fd(longueur)
            pd()


############################ Exécution du TP #################################

## On crée d'abord le brin complémentaire grâce à la fonction antisens puis
## on colle les deux brins afin de nous apprêter à la recherche de gènes.
## On appelle le tout "ADN".
           
ADN = adn + antisens(adn)

## Deuxièmement, nous allons parcourir "ADN" pour trouver tous les codons de
## départ et tous les codons de terminaison grâce à nos fonctions trouveDebut
## et trouveFin. On stocke le tout dans des variables pour faciliter la
## compréhension.

ADNPositionsDebut = trouveDebut(ADN) 
ADNPositionsFin = trouveFin(ADN) 

## Ensuite, nous allons nous baser sur ces deux derniers tableaux pour générer
## le tableau de tuples contenant la liste des gènes trouvés sur notre "ADN".
## Encore une fois, nous stockons le tout dans une variable pour faciliter la
## compréhension.

ADNPositionsGenes = trouveGene(ADNPositionsDebut,ADNPositionsFin) 

## Puis, grâce à notre fonction intermédiaire lesGenes, nous allons générer un 
## tableau contenant les sous-chaines de caractères du brin "ADN".

ADNGenes = lesGenes([(38, 314), (74, 314), (402, 723), (751, 811)],ADN)  

## Et nous sélectionnons nos gènes que nous stockons dans de nouvelles varia-
## bles temporaires, pour ensuite les transformer en ARN grâce à notre fonction
## transcrire.

temp1 = ADNGenes[0]
temp2 = ADNGenes[1]
temp3 = ADNGenes[2]
temp4 = ADNGenes[3]

gene1 = transcrire(temp1)
gene2 = transcrire(temp2)
gene3 = transcrire(temp3)
gene4 = transcrire(temp4)


## Grâce à notre fonction fromArnToProteine, nous allons afficher chaque 
## gène (sous forme d'ARN) en protéine à l’aide du premier tableau associatif
## fourni dans le fichier de départ.

print(fromArnToProteine(gene1))
print(fromArnToProteine(gene2))
print(fromArnToProteine(gene3))
print(fromArnToProteine(gene4))
    

## Enfin, grâce à nos fonctions lettresDesCodons et carre, nous allons
## afficher les protéines à l’aide de la tortue en utilisant le deuxième 
## tableau associatif, chaque acide aminé étant représenté par une lettre.

## On se met dans la bonne position et on représente les protéines du premier
## gène.
goto(-160,250)
lettresDesCodons(gene1)

## On fait de même pour le deuxième gène.

pu()
fd(30)
lt(90)
bk(20)
pd()

lettresDesCodons(gene2)

## Même chose pour le troisième.

pu()
fd(30)
lt(90)
bk(80)
pd()

lettresDesCodons(gene3)

## Et enfin, même chose pour le quatrième.

pu()
fd(30)
lt(90)
bk(20)
pd()

lettresDesCodons(gene4)
  
## Pour terminer, on va tracer nos carrés grâce à la fonction carre.

# D'abord pour le premier gène.

pu()
lt(90)
goto(-170,260)
pd()

carre(20,91)

# Pour replacer la tortue, comme le mécanisme est identique à chaque fois,
# nous allons définir une fonction "replacer" à cet effet.

def replacer(nombre):
    pu()
    rt(90)
    fd(30)
    lt(90)
    bk(nombre)
    pd()

# Ensuite, le deuxième.

replacer(20)
carre(20,79)
    
    
# Puis, le troisième.

replacer(80)
carre(20,106)

# Et enfin, le quatrième.

replacer(20)
carre(20,19)

################################## FIN #######################################

replacer(20)
