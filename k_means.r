# Norme euclidienne
norme <- function(x){ return(sqrt(sum(x^2))) }

#################
#### K-Means ####
#################

Kmeans <- function(data, K, npas_max=50, afficher=F){
  
  # Initialisation
  I_centres = sample.int(nrow(data),size = K)
  centres = (data[I_centres,])
  
  ## Repetition de l'algorithme jusqu'au critere d'arret : 
  # - Soit le nombre de pas maximal est atteint, 
  # - soit les clusters obtenus ne changent plus d'une iteration a l'autre et alors diff = 1
  diff = 0
  npas = 0
  
  while(npas < npas_max && diff != 1){
    npas = npas + 1
    
    # Partition en K clusters
    nouveaux_clusters_choisis = as.factor(apply(data, 1, function(point){    # Pour chaque point,
      # On calcule la distance a chaque centre 
      distance_aux_centres = apply(centres,1, function(x){ return(norme(x-point))  })
      # Puis on prend la plus petite (avec which.min)
      return(as.numeric(which.min(distance_aux_centres)))     
    }))
    
    # On affiche les resultats si demandes
    if(afficher==T){
      graphique= ggplot(data, aes(x=X1,y=X2))+
        geom_point(color=nouveaux_clusters_choisis)+
        geom_point(data=centres, aes(x=X1, y=X2), color=as.factor(1:4), size=5)+
        ggtitle(paste("Algorithme K-means - Pas",npas))
      print(graphique)
    }
    
    # On verifie si les clusters obtenus ont change ou non (sauf au premier pas)
    if(npas == 1){clusters_choisis = nouveaux_clusters_choisis}
    else {if(all(clusters_choisis == nouveaux_clusters_choisis)){diff = 1}}
    clusters_choisis = nouveaux_clusters_choisis

    
    # Choix des nouveaux centres des clusters
    centres = aggregate(data,list(clusters_choisis), mean)[,2:3]
    
  }
  
  return(clusters_choisis)
  
}

npas_Kmeans <- function(data, K, npas_max=50, n=10){
  liste_npas = c()
  
  for(i in 1:n){
    print(paste("simulation",i))
    # Initialisation
    I_centres = sample.int(nrow(data),size = K)
    centres = (data[I_centres,])
    
    ## Repetition de l'algorithme jusqu'au critere d'arret : 
    # - Soit le nombre de pas maximal est atteint, 
    # - soit les clusters obtenus ne changent plus d'une iteration a l'autre et alors diff = 1
    npas = 0          
    diff = 0
    while(npas < npas_max && diff != 1){
      npas = npas + 1
      
      # Partition en K clusters
      nouveaux_clusters_choisis = as.factor(apply(data, 1, function(point){    # Pour chaque point,
        # On calcule la distance a chaque centre 
        distance_aux_centres = apply(centres,1, function(x){ return(norme(x-point))  })
        # Puis on prend la plus petite (avec which.min)
        return(as.numeric(which.min(distance_aux_centres)))     
      }))
      
      
      # On verifie si les clusters obtenus ont change ou non (sauf au premier pas)
      if(npas == 1){clusters_choisis = nouveaux_clusters_choisis}
      else {if(all(clusters_choisis == nouveaux_clusters_choisis)){diff = 1}}
      clusters_choisis = nouveaux_clusters_choisis
      
      
      # Choix des nouveaux centres des clusters
      centres = aggregate(data,list(clusters_choisis), mean)[,2:3]
      
    }
    liste_npas = c(liste_npas, npas)
  }
  ggplot(data.frame(x=liste_npas), aes(x))+geom_histogram(aes(y=..count../sum(..count..)), binwidth=2)+
    ggtitle(paste("Histogramme du nombre de pas de K-means sur",n,"essais"))+
    xlab("Nombre de pas avant convergence")+
    ylab("Proportion")
  
}

###############
#### Bonus ####
###############

# Calcul de l'inertie (somme des distances de chaque point au centre du cluster) 

# S = sum(apply(data - centres[clusters_choisis,], 1, norme))
# S 

