Paradevs-learning
=================

Transformation du graphe de modèles PDEVS en graphe pondérés.

On considère que chaque modèle atomique possède un unique port d'entrée et un
unique port de sortie. Les événements ne transportent pas de données. Le graphe
se compose de N types de modèles où N est inférieur au nombre global de
modèles du graphe.

Etape 1 : (-D) - phase de génération de séquence d'apprentissage
 * chaque type de modèle fait l'objet d'une génération de séquences ;
 * on connecte à l'entrée du modèle un générateur d'événements aléatoires ;
 * une séquence est une chaîne de caractères composée "a|O|9" et "a||INF" où la
   première lettre est l'identifiant du modèle, la lettre O si c'est un
   événement de sortie (sinon rien) puis la durée écoulée depuis le
   précédent événement de sortie ou INF si on est dans un état à durée infinie
 * la séquence générée est fonction de la durée de simulation
 * la séquence est envoyée sur la console

Etape 2 : (-L x n i) - construction des matrices initiales (emit et trans) et
      apprentissage
 * une première phase consiste à construire les matrices emit et trans des
   HMM en fonction des symboles de la séquence et du nombre d'états cachés (n)
 * la séquence doit être contenu dans un fichier nommé "seq"
 * un seul type de modèle est traité (x)
 * puis la séquence initiale est découpée en sous-séquences de longueur 20 dans
   un fichier seq-xxx.input afin de constituer un ensemble d'apprentissage ;
   chaque sous-séquence se termine par un symbole de fin (END)
 * une séquence est composée d'observations ; une observation est un couple
   (durée|type) où type = 1 si événement de sortie, type = 2 si état infini ou
   type = 3 si symbole de fin de séquence
 * puis la deuxième phase consiste en l'apprentissage proprement dit ; le
   nombre d'itérations est fixé (i)
 * les fichiers seq-xxx-init.yyy contiennent les probabilités initiales
 * les fichiers seq-xxx-result.yyy contiennent les probabilités après
   l'apprentissage

Etape 3 : (-G n ...) - génération des poids
 * chaque modèle est remplacé par son HMM équivalent
 * on simule les modèles sur la période [0, n]
 * à l'aide des matrices et symboles, on calcule les dates d'envoi des
   événements
 * deux cas de figure :
   - l'état est de type "événement de sortie" (OUT) : la date courante et la
     date d'envoi sont égales à la dernière date courante + la durée de l'état
   - l'état est infini : la date courante est égale à la date d'arrivée
     du prochain événement d'entrée (pas d'envoi d'événement de sortie)
 * on calcule la moyenne des durées entre 2 événements ; ce qui constituera
   le poids du modèle

Etape 4 : (-M x) - calcul de la moyenne à partir de la séquence initiale
 * à partir de la séquence initiale et en précisant le nom du modèle (x),
   la moyenne des durées entre 2 événements est calculée
 * cela permet de voir l'écart avec l'apprentissage

===============================================================================

Exemple :
A, B => D
B, C => E
D, E => F
E -> G
F, G => H

Tous les modèles possèdent la même dynamique, seul le graphe de connexion
change.

Plusieurs stratégies d'apprentissage :
* Si on simule la totalité du graphe pour générer la séquence initiale, avec 10
états cachés, le HMM est très proche.
* Si on simule seulement un modèle avec un générateur d'événements aléatoires et
que l'on utilise le HMM appris sur cette séquence, tous les modèles possèdent
la même moyenne.
* Si on simule un modèle sans entrée pour les modèles A,B et C + on simule un
modèle à deux entrées aléatoires pour D, E, F et H + on simule un modèle à une
entrée aléatoires pour G, on obtient un meilleur résultat mais divergent pour
les modèles avec entrées
* Si on simule un modèle sans entrée pour les modèles A, B et C, on apprends et
on remplace A, B et C et on simule A, B, C et D (où D est le modèle DEVS et E
est identique à D, pas besoin de le simuler)
; et on itére le processus => je pense que l'on approche d'une bonne solution
sans tout simuler ==> A TESTER !