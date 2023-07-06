
Algorithme : read
Entrée : nom_fichier, timestep, N_atomes, N_conf, pos, labels, selection
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
	
	// Position de départ dans le fichier
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
	
	
	allouer(pos, N_conf * N_atomes * 3)
	
	
	// READING
	aller_fichier(fichier, debut_fichier)
	Pour c de 0 à N_conf - 1
	Faire
		Pour a de 0 à N_atomes - 1
		Faire
			Si selection[a]
			Alors
				lire_format(fichier, "%d %*s ", indice)
				lire_format(fichier, "%f %f %f %*f %*f %*f %*f", pos[c][indice][0], pos[c][indice][1], pos[c][indice][2])
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
Entrée : nom_fichier, timestep, N_conf, steps, bounds, N_selection, indices, pos, charges
Description :
	nom_fichier est utilisé pour ouvrir et lire le fichier des configurations
	Le timestep correspond au timestep à partir duquel on doit commencer à lire les configurations
	N_conf est le nombre de configurations à traiter
	steps est le tableau qui contiendra les numéros des configurations (en itérations de simulation)
	bounds est le tableau qui contiendra les limites de la boîte de simulation
	N_atomes est le tableau qui contiendra les nombres d'atomes pour chaque configuration
	N_selection est le nombre d'atomes correspondant à la séléction
	pos est le tableau qui contiendra les positions de chacun des atomes pour chaque configuration
	charges est le tableau qui contiendra les charges de chacun des atomes pour chaque configuration
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
			# Note d'implémentation : ajuster la taille du tableau steps
			steps[N_conf - 1] <- step
		Sinon Si texte = "ITEM: NUMBER OF ATOMS"
		Alors
			lire_format(fichier, "%d\n", N_atomes)
		Sinon Si texte = "ITEM: BOX BOUNDS pp pp pp\n"
		Alors
			# Note d'implémentation : ajuster la taille du tableau bounds
			Pour d de 0 à 2
			Faire
				lire_format(fichier, "%lf %lf\n", bounds[N_conf - 1][2 * d], bounds[N_conf - 1][2 * d + 1])
			Fin Pour
		Sinon Si texte = "ITEM: ATOMS id element x y z ... q"
		Alors
			# Note d'implémentation : ajuster la taille des tableaux N_selection, indices, pos et charges
			N_selection[N_conf - 1] <- 0
			Pour a de 0 à N_atomes - 1
			Faire
				lire_format(fichier, "%d %s", indice, element)
				Si element != "C"
				Alors
					lire_ligne(fichier)
					passer
				Fin Si
				N_selection[N_conf - 1] <- N_selection[N_conf - 1] + 1
				indices[N_conf - 1][N_selection[N_conf - 1] - 1] <- indice
				lire_format(fichier, "%lf %lf %lf %*f %*f %*f %lf\n", pos[N_conf - 1][N_selection[N_conf - 1] - 1][0:3], charges[N_conf - 1][N_selection[N_conf - 1] - 1])
			Fin Pour
			# Note d'implémentation : ajuster la taille des tableaux modifiés
		Sinon
			erreur("Erreur de lecture")
		Fin Si
	Fin Tant que
Fin


Algorithme : bonds
Données : pos, N_conf, N_atomes, selection
Entrée : N_bonds, bonds, R, D_xy, D_z
Début
	allouer(N_bonds, N_conf * N_atomes)
	allouer(bonds, N_conf * N_atomes)
	Pour c de 0 à N_conf
	Faire
		Pour i de 0 à N_atomes
		Faire
			Si selection[i]
			Alors
				Pour j de i + 1 à N_atomes
				Faire
					Si selection[j]
					Alors
						Si - D_z <= pos[c][j][2] - pos[c][i][2] <= D_z et
						   - D_xy <= (pos[c][j][0] - pos[c][i][0])² + (pos[c][j][1] - pos[c][i][1])² - R² <= D_xy
						Alors
							N_bonds[c][i] <- N_bonds[c][i] + 1
							reallouer(bonds[c][i], N_bonds[c][i])
							bonds[c][i][N_bonds - 1] <- j + 1
							
							N_bonds[c][j] <- N_bonds[c][j]
							reallouer(bonds[c][j], N_bonds[c][j])
							bonds[c][j][N_bonds - 1] <- i + 1
						Fin Si
					Fin Si
				Fin Pour
			Fin Si
		Fin Pour
	Fin Pour
Fin


Algorithme : bonds_2
Données : N_conf, N_selection, bounds, indices, pos, charges
Entrée : N_bonds, bonds, R, delta
Début
	# Note d'implémentation : allouer N_bonds, bonds
	# Note d'implémentation : bonds est alloué avec 4 emplacements par atome (le nombre maximal théorique pour le carbone)
	Pour c de 0 à N_conf - 1
	Faire
		Pour i de 0 à N_selection[c] - 1
		Faire
			Pour j de i + 1 à N_selection[c] - 1
			Faire
				r2 <- 0
				Pour d de 0 à 3 - 1
				Faire
					length <- (bounds[c][2 * d + 1] - bounds[c][2 * d])
					diff <- pos[c][j][d] - pos[c][i][d]
					Si diff < - length / 2.
					Alors
						diff <- diff + length
					Sinon Si length / 2. < diff
					Alors
						diff <- diff - length
					Fin Si
					r2 <- r2 + diff * diff
				Fin Pour
				
				Si (R - delta) * (R - delta) <= r2 <= (R + delta) * (R + delta)
				Alors
					N_bonds[c][i] <- N_bonds[c][i] + 1
					bonds[c][i][N_bonds - 1] <- j
					
					N_bonds[c][j] <- N_bonds[c][j]
					bonds[c][j][N_bonds - 1] <- i
				Fin Si
			Fin Pour
		Fin Pour
		# Note d'implémentation : ajuster la taille des emplacements de bonds
	Fin Pour
Fin


Algorithme : layers
Données : pos, N_conf, N_atomes, selection
Entrée : layers, H, D_z
Début
	allouer(layers, N_conf * N_atomes)
	Pour c de 0 à N_conf
	Faire
		Pour i de 0 à N_atomes
		Faire
			haut <- faux
			bas <- faux
			layers[c][i] <- faux
			Si selection[i]
			Alors
				Pour j de i + 1 à N_atomes
				Faire
					Si selection[j]
					Alors
						Si non haut et - D_z <= pos[c][j][2] - pos[c][i][2] - H <= D_z
						Alors
							haut <- vrai
						Fin Si
						Si non bas et - D_z <= pos[c][j][2] - pos[c][i][2] + H <= D_z
						Alors
							bas <- vrai
						Fin Si
						Si haut et bas
						Alors
							layers[c][i] <- vrai
							stop
						Fin Si
					Fin Si
				Fin Pour
			Fin Si
		Fin Pour
	Fin Pour
Fin


Algorithme : electrodes
Données : N_conf, N_selection, pos
Entrée : electrodes, limit
Début
	# Note d'implémentation : allouer le tableau electrodes
	Pour c de 0 à N_conf - 1
	Faire
		Pour a de 0 à N_selection[c] - 1
		Faire
			electrodes[c][a] <- pos[c][a][2] >= limit
		Fin Pour
	Fin Pour
Fin


Algorithme : histograms
Données : N_conf, N_selection, layers, electrodes, charges
Entrée : histograms, N_types
Début
	# Note d'implémentation : allouer les tableaux histograms et N_types
	Pour c de 0 à N_conf - 1
	Faire
		Pour a de 0 à N_selection[c] - 1
		Faire
			Si layers[c][a] = 0 et electrodes[c][a] = 0
			Alors
				t <- 0
			Sinon Si layers[c][a] = 0 et electrodes[c][a] = 1
			Alors
				t <- 1
			Sinon Si layers[c][a] = 1 et electrodes[c][a] = 0
			Alors
				t <- 2
			Sinon
				t <- 3
			Fin Si
			histograms[c][t] <- histograms[c][t] + charges[c][a]
			N_types[c][t] <- N_types[c][t] + 1
		Fin Pour

		Pour t de 0 à 3 - 1
		Faire
			histograms[c][t] <- histograms[c][t] / N_types[c][t]
		Fin Pour
	Fin Pour
Fin
