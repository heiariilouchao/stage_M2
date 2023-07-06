
Algorithme : read
Entrée : nom_fichier, timestep, N_atomes, N_conf, upos, labels, selection
Début
	fichier <- ouvrir_fichier(nom_fichier)
	
	
	// SKIPPING
	texte <- ""
	step <- 0
	Tant que step != timestep
	Faire
		Tant que texte != "ITEM: TIMESTEP\n"
		Faire
			texte <- lire_ligne(fichier)
		Fin Tant que
		step <- lire_entier(fichier)
	Fin Tant que
	
	
	// PRE-READING
	// Nombre d'atomes
	Tant que texte != "ITEM: NUMBER OF ATOMS\n"
	Faire
		texte <- lire_ligne(fichier)
	Fin Tant que
	N_atomes <- lire_entier(fichier)
	
	// Position de départ
	Tant que texte != "ITEM: ATOMS ...\n"
	Faire
		texte <- lire_ligne(fichier)
	Fin Tant que
	debut_fichier <- pos_fichier(fichier)

	// Sélection des atomes
	Pour a de 0 à N_atomes - 1
	Faire
		lire_format(fichier, "%d %s", indice, label)
		Pour l de 0 à N_labels - 1
		Faire
			Si label = labels[l]
			Alors
				selection[indice] <- vrai
				stop
			Fin Si
		Fin Pour
	Fin Pour

	// Nombre de configurations
	aller_fichier(fichier, debut_fichier)
	N_conf <- 1
	Tant que non fin_fichier(fichier)
	Faire
		Si lire_ligne(fichier) = "ITEM: ATOMS ...\n"
		Alors
			N_conf <- N_conf + 1
		Fin Si
	Fin Tant que
	
	
	allouer(upos, N_conf * N_atomes * 3)
	
	
	// READING
	aller_fichier(fichier, debut_fichier)
	Pour c de 0 à N_conf - 1
	Faire
		Pour a de 0 à N_atomes - 1
		Faire
			Si selection[a]
			Alors
				lire_format(fichier, "%d %*s ", indice)
				lire_format(fichier, "%*f %*f %*f %lf %lf %lf %*f", upos[c][indice][0], upos[c][indice][1], upos[c][indice][2])
			Fin Si
		Fin Pour
		Si c < N_conf - 1
		Alors
			Faire
				texte <- lire_ligne(fichier)
			Tant que texte != "ITEM: ATOMS ...\n"
		Fin Si
	Fin Pour
Fin


Algorithme : read_2
Entrée : nom_fichier, timestep, N_conf, N_selection, selection, indices, upos
Description :
	nom_fichier est utilisé pour ouvrir et lire le fichier des configurations
	timestep correspond au timestep à partir duquel on doit commencer à lire les configurations
	N_conf est le nombre de configurations à traiter
	N_selection est le nombre d'atomes correspondant à la séléction
	selection est une chaine de caractères correspondant à la liste des éléments à sélectionner
	indices est le tableau qui contiendra les indices des atomes
	upos est le tableau qui contiendra les positions déroulées de chacun des atomes pour chaque configuration
Début
	fichier <- ouvrir_fichier(nom_fichier)
	
	N_conf <- 0
	Tant que non fin_fichier(fichier)
	Faire
		texte <- lire_ligne(fichier)
		Si texte = "ITEM: TIMESTEP\n"
		Alors
			lire_format(fichier, "%d\n", step)
			Tant que step < timestep
			Faire
				Faire
					texte <- lire_ligne(fichier)
				Tant que texte != "ITEM: TIMESTEP\n"
				Fin Tant que
				lire_format(fichier, "%d\n", step)
			Fin Tant que
			N_conf <- N_conf + 1
		Sinon Si texte = "ITEM: NUMBER OF ATOMS"
		Alors
			lire_format(fichier, "%d\n", N_atomes)
		Sinon Si texte = "ITEM: BOX BOUNDS pp pp pp\n"
		Alors
			Pour d de 1 à 3
			Faire
				lire_ligne(fichier)
			Fin Pour
		Sinon Si texte = "ITEM: ATOMS id element x y z ... q"
		Alors
			# Note d'implémentation : ajuster la taille des tableaux N_selection, indices et upos
			N_selection[N_conf - 1] <- 0
			Pour a de 0 à N_atomes - 1
			Faire
				lire_format(fichier, "%d %s", indice, element)
				Si non trouver_sous_chaine(selection, element)
				Alors
					lire_ligne(fichier)
					passer
				Fin Si
				N_selection[N_conf - 1] <- N_selection[N_conf - 1] + 1
				indices[N_conf - 1][N_selection[N_conf - 1] - 1] <- indice
				lire_format(fichier, "%*f %*f %*f %lf %lf %lf %*f\n", upos[N_conf - 1][N_selection[N_conf - 1] - 1][0:3])
			Fin Pour
			# Note d'implémentation : ajuster la taille des tableaux modifiés
		Sinon
			erreur("Erreur de lecture")
		Fin Si
	Fin Tant que
Fin


Algorithme : msd
Données : upos, N_conf, N_atomes, selection
Entrée : msd
Début
	Pour tau de 1 à N_conf - 1
	Faire
		Pour c de 0 à N_conf - 1 - tau
		Faire
			Pour a de 0 à N_atomes - 1
			Faire
				Si selection[a]
				Alors
					dist <- 0
					Pour d de 0 à 2
					Faire
						dist <- dist + carre(upos[c + tau][a][d] - upos[c][a][d])
					Fin Pour
					msd[tau] <- msd[tau] + racine(dist)
				Fin Si
			Fin Pour
		Fin Pour
		msd[tau] <- 1 / (N_conf - 1 - tau) * 1 / N_atomes * msd[tau]
	Fin Pour
Fin


Algorithme : msd_2
Données : N_conf, N_selection, indices, upos
Entrée : msd
Début
	Pour tau de 1 à N_conf - 1
	Faire
		Pour c de 0 à N_conf - 1 - tau
		Faire
			somme <- 0.
			Pour a de 0 à N_selection[c] - 1
			Faire
				dist <- 0.
					Pour d de 0 à 3 - 1
					Faire
						dist <- dist + carre(upos[c + tau][a][d] - upos[c][a][d])
					Fin Pour
				somme <- somme + racine(dist)
			Fin Pour
			msd[tau] <- msd[tau] + somme / N_selection[c]
		Fin Pour
		msd[tau] <- msd[tau] / (N_conf - 1 - tau)
	Fin Pour
Fin
