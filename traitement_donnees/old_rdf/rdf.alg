Algorithme : parse_elements
Données : sep
Entrée : elements, N_elements, labels
Début
	Si elements[0] = sep ou elements[taille_chaine(elements) - 1] = sep
	Alors
		erreur
	Fin Si
	N_elements <- 1
	debut <- 0
	Pour c de 0 à taille_chaine(elements) - 1
	Faire
		Si elements[c] = sep
		Alors
			# Note d'implémentation : allouer labels[N_elements - 1] avec (c - debut) + 1 caractères
			copier_n_chaine(labels + N_elements - 1, elements + debut, (c - debut))
			N_elements <- N_elements + 1
			debut <- c + 1
		Fin Si
	Fin Pour
	# Note d'implémentation : allouer labels[N_elements - 1] avec (taille_chaine(elements) - debut) + 1 caractères
	copier_n_chaine(labels + N_elements - 1, elements + debut, (taille_chaine(elements) - debut))
Fin



Algorithme : read
Entrée : nom_fichier, timestep, N_conf, bounds, labels, N_elements, elements, N_selection, indices, pos
Description :
	nom_fichier est utilisé pour ouvrir et lire le fichier des configurations
	Le timestep correspond au timestep à partir duquel on doit commencer à lire les configurations
	N_conf est le nombre de configurations à traiter
	bounds est le tableau qui contiendra les limites de la boîte de simulation
	N_selection est le nombre d'atomes correspondant à chaque type pour chaque configuration
	indices est le tableau qui contiendra les indices des atomes sélectionnés pour chaque configuration
	pos est le tableau qui contiendra les positions de chaque atome pour chaque configuration
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
			# Note d'implémentation : ajuster la taille du tableau bounds
			Pour d de 0 à 3 - 1
			Faire
				lire_format(fichier, "%lf %lf\n", bounds[N_conf - 1][2 * d], bounds[N_conf - 1][2 * d + 1])
			Fin Pour
		Sinon Si texte = "ITEM: ATOMS id element..."
		Alors
			# Note d'implémentation : ajuster la taille des tableaux N_selection et pos
			# Note d'implémentation : allouer et initialiser à zéro les entrées de N_selection
			Pour a de 0 à N_atomes - 1
			Faire
				lire_format(fichier, "%d %s", indice, element)
				Si non chaine_dans_chaine(labels, element)
				Alors
					lire_ligne(fichier)
					passer
				Fin Si
				Pour e de 0 à N_elements - 1
				Faire
					Si element = elements[e]
					Alors
						type <- e
					Fin Si
				Fin Pour
				N_selection[N_conf - 1][type] <- N_selection[N_conf - 1][type] + 1
				indices[N_conf - 1][type][N_selection[N_conf - 1][type] - 1] <- indice
				Si texte = "ITEM: ATOMS id element x y z xu yu zu q"
				Alors
					lire_format(fichier, "%lf %lf %lf %*f %*f %*f %*f\n", pos[N_conf - 1][type][N_selection[N_conf - 1][type] - 1][0:3])
				Sinon Si texte = "ITEM: ATOMS id element x y z xu yu zu"
				Alors
					lire_format(fichier, "%lf %lf %lf %*f %*f %*f\n", pos[N_conf - 1][type][N_selection[N_conf - 1][type] - 1][0:3])
				Fin Si
			Fin Pour
			# Note d'implémentation : ajuster la taille des tableaux modifiés
		Sinon
			erreur("Erreur de lecture")
		Fin Si
	Fin Tant que
Fin


Algorithme : transpose
Données : N_elements
Entrée : e1, e2
Sortie : i
Début
	somme <- 0
	Pour e de 0 à e1 - 1
	Faire
		somme <- somme + (N_elements - e)
	Fin Pour
	Retourner somme + e2 - e1
Fin


Algorithme : rdf
Données : N_conf, N_elements, N_selection, bounds, pos
Entrée : N_bins, hist, r, RDF
Début
	N_pairs <- N_elements * (N_elements + 1) / 2
	# Note d'implémentation : allouer hist, r et RDF avec N_bins
	Pour b de 0 à N_bins
	Faire
		r[b] <- b * delta
	Fin Pour

	Pour c de 0 à N_conf - 1
	Faire
		Pour e1 de 0 à N_elements - 1
		Faire
			Pour e2 de e1 à N_elements - 1
			Faire
				Pour a1 de 0 à N_selection[c][e1] - 1
				Faire
					Pour a2 de a1 à N_selection[c][e2] - 1
					Faire
						r2 <- 0.
						Pour d de 0 à 3 - 1
						Faire
							length <- (bounds[c][2 * d + 1] - bounds[c][2 * d])
							diff <- pos[c][e2][a2][d] - pos[c][e1][a1][d]
							Si diff < - length / 2.
							Alors
								diff <- diff + length
							Sinon Si length / 2. < diff
							Alors
								diff <- diff - length
							Fin Si
							r2 <- r2 + diff * diff
						Fin Pour
						Si r2 / delta < N_bins
						Alors
							somme[r2 / delta] <- somme[r2 / delta] + 2
						Fin Si
					Fin Pour
				Fin Pour

				N <- N_selection[c][e1]
				Si e1 = e2
				Alors
					N <- N * (N - 1)
				Sinon
					N <- N * N_selection[c][e2]
				Fin Si

				Pour b de 0 à N_bins - 1
				Faire
					hist[transpose(e1, e2)][b] <- somme[b] / N
			Fin Pour
		Fin Pour
	Fin Pour

	const <- 1. / (4. / 3. * PI / V * N_conf)
	Pour p de 0 à N_pairs - 1
	Faire
		Pour b de 0 à N_bins - 1
		Faire
			RDF[p][b] <- const * hist[p][b] / ((r[b] + delta)^3 - r[b]^3)
			r[b] <- r[b] + delta / 2.
		Fin Pour
	Fin Pour
Fin




0 1 2 3
00 01 02 03
11 12 13
22 23
33
-> i(e1, e2) = somme_{0}^{e1 - 1}(N_elements - j) + (e2 - e1)
