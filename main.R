# Chemin du dossier ou sont enregistres les 2 fichiers .R et la base "EMGaussian_csv" 

setwd("")                ##### A MODIFIER

# Importations des librairies specifiques

require(mixtools)
require(ggplot2)
source("k_means.r")


#################
#### Donnees ####
#################

data = read.csv("EMGaussian.csv", sep=" ",header=F)
colnames(data) <- c("X1","X2")

####################
#### Parametres ####
####################

K = 4
n = nrow(data)

#########################################################
#### Initialisation de l'algorithme EM avec K-means #####
#########################################################

initialisation <- function(data, K, npas_max, afficher = F){
  
  clusters = Kmeans(data, K, npas_max, F)

  # Moyennes des clusters
  moyennes = aggregate(data,list(clusters), mean)[,2:3] 
  # La fonction aggregate calcule la moyenne des clusters (c'est a dire les coordonnees(X1,X2) des barycentres)
  # selon le cluster (1,2,3 ou 4)
  
  # Matrice de variance-covariance des clusters
  length_clusters = c(table(clusters))
  variances = t(apply((data - moyennes[clusters,])/sqrt(length_clusters[clusters]), 1, function(x) as.matrix(x)%*%t(x) ))
  matrices_variances = aggregate(variances, list(clusters), sum)[,-1]
  
  # Variances scalaires dans le cas des matrices proportionnelles (on calcule la variance empirique sur chaque cluster : somme(t(x_i - mu_i)*(x_i - mu_i) )
  variances_scalaires = apply((data - moyennes[clusters,])/sqrt(length_clusters[clusters]), 1, function(x) t(x)%*%as.matrix(x) )
  liste_variances_scalaires = aggregate(variances_scalaires, list(clusters), sum)[,-1]
  
  # Calcul des poids alpha_k
  
  alpha = c(table(clusters)/n)
  
  # Calcul de la log-vraisemblance du modele a donnees manquantes
  
  log_vraisemblance = 0
  s = 0
  
  for(i in 1:n){
    for(l in 1:K){
      x_mu_l = t(as.matrix(data[i,] - moyennes[l,]))
      sigmal_t =  matrix(as.numeric(matrices_variances[l,]),2,2)
      s = s + (alpha[l]/( (2*pi)*sqrt(det(sigmal_t)) )) * exp(-0.5 *t(x_mu_l)%*%solve(sigmal_t)%*%x_mu_l) 
    }
    log_vraisemblance = log_vraisemblance + log(s)
  }
  
  if (afficher == T){
    graphique= ggplot(data, aes(x=X1,y=X2,color=clusters))+
      geom_point()+
      geom_point(data=moyennes, aes(x=X1, y=X2), size=5, color="gray30")+
      stat_ellipse(type="norm", level=0.80)+
      ggtitle(paste("Resultat K-means pour",K,"clusters"))
    print(graphique)
  }
  
  return(list("log_vraisemblance" = log_vraisemblance, "alpha" = alpha, "moyennes" = moyennes, "matrices_variances" = matrices_variances , "liste_variances_scalaires" = liste_variances_scalaires))
}

########################   
#### Algorithme EM #####
########################   

Algo_EM <- function(data, K, npas_max, maxiter, epsilon, afficher=T, afficher_ellipse=F, afficher_cercle=F){
  
  P = matrix(nrow = n, ncol = K)
  delta = 1 # valeur arbitraire superieure a epsilon
  
  # Initialisation : a chaque fois que l'algorithme EM est lance, une nouvelle initalisation a lieu
  
  init = initialisation(data, K, npas_max, F)
    
  log_vraisemblance = init$log_vraisemblance
  moyennes = init$moyennes
  matrices_variances = init$matrices_variances
  liste_variances_scalaires = init$liste_variances_scalaires
  alpha = init$alpha
  
  t = 0
  
  # Double critere d'arret
  
  while (t < maxiter & delta > epsilon) {
    t = t+1
    ## Etape maximisation de l'algorithme
    
    # Calcul des p_ik(t), i.e les probabilites que une observation X_i appartienne a la classe k
    
    new_log_vraisemblance = 0 # Dans le meme temps on calcule la nouvelle log-vraisemblance du modele a donnees manquantes
    
    for(i in 1:n){
      
      #denominateur 
      den = 0
      for(l in 1:K){
        x_mu_l = t(as.matrix(data[i,] - moyennes[l,]))
        sigmal_t =  matrix(as.numeric(matrices_variances[l,]),2,2)
        den = den + (alpha[l]/( (2*pi)*sqrt(det(sigmal_t)) )) * exp(-0.5 *t(x_mu_l)%*%solve(sigmal_t)%*%x_mu_l) 
      }
      new_log_vraisemblance = new_log_vraisemblance + log(den)
      
      # p_ik(t)
      for(k in 1:K){
        x_mu_k=t(as.matrix(data[i,] - moyennes[k,]))
        sigmak_t =  matrix(as.numeric(matrices_variances[k,]),2,2)
        #numerateur
        num = (alpha[k]/( (2*pi)*sqrt(det(sigmak_t)) )) * exp(-0.5 *t(x_mu_k)%*%solve(sigmak_t)%*%x_mu_k  ) 
        P[i, k] = num / den
      }
    }
    
    delta = abs(new_log_vraisemblance - log_vraisemblance)
    log_vraisemblance <- new_log_vraisemblance
    
    # alpha_k(t+1)
    alpha = apply(P,2,mean)
    
    # mu_k(t+1)
    for(k in 1:K){
      moyennes[k,]  = c(0,0)
      for(i in 1:n){
        moyennes[k,] = moyennes[k,] + (P[i,k] / sum(P[,k])) * data[i,]
      } 
    }
    
    #sigma_k(t+1)
    for(k in 1:K){
      matrices_variances[k,] = c(0,0,0,0)
      for(i in 1:n){
        
        x_mu_k=t(as.matrix(data[i,] - moyennes[k,]))
        
        matrices_variances[k,] = matrices_variances[k,] + c(t( (P[i,k] / sum(P[,k])) * x_mu_k%*%t(x_mu_k) ))
        
      }   
    }
    
    #sigma(t+1) scalaire
    liste_variances_scalaires = c(0,0,0,0)
    for(k in 1:K){
      
      for(i in 1:n){
        
        x_mu_k=t(as.matrix(data[i,] - moyennes[k,]))
        
        liste_variances_scalaires[k] = liste_variances_scalaires[k] + (P[i,k] / (2*sum(P[,k]))) * t(x_mu_k)%*%x_mu_k 
        
      }   
    }
    
    if(afficher==T){
      # La ligne ci-dessous correspond a la regle du Maximum a Posteriori : on assigne chaque observation au cluster dont elle a la plus forte probabilit? d'appartenance
      clusters = factor(apply(P, 1, which.max))
      
      if(afficher_cercle==T){
        ## Creation des zones de confiance a 80% pour la matrice de variance scalaire (zone de confiance = disque)
        df_ell <- data.frame()
        for(cluster in 1:4){
          # Pour chaque cluster, on recupere sa moyenne et sa matrice de variance covariance
          mu = c(as.numeric(moyennes[cluster,]))
          sigma_scalaire = liste_variances_scalaires[cluster]*diag(2)
          # On Calcule les coordonnees des points de chaque cercle de confiance grace a la fonction ellipse du package mixtools
          df_ell <- rbind(df_ell, cbind(as.data.frame(mixtools::ellipse(mu, sigma_scalaire, npoints=250, alpha=0.2, draw=F)), clusters=as.factor(cluster) ))
        }
        colnames(df_ell) <- c("x","y","clusters")
      }
      
      # On affiche les zones de confiance (ellipses et cercles) correspondant aux cas scalaire et cas general si demande (avec affiche_ellipse et affiche_cercle)
      graphique= ggplot(data, aes(x=X1,y=X2,color=clusters))+
        geom_point()+
        geom_point(data=moyennes, aes(x=X1, y=X2), size=5, color="gray30")+
        ggtitle(paste("Algorithme EM - etape", t)) 
      
      if(afficher_ellipse==T){
        graphique = graphique+stat_ellipse(type="norm", level=0.80)
      }
      if(afficher_cercle==T){
        graphique = graphique+geom_path(data=df_ell, aes(x,y,color=clusters), size=1, linetype=2)
      }
      
      
      print(graphique)
    }
    
    print(paste("etape",t,":","delta = ", delta ))
    
  }
  
}


#################################
#### Execution des fonctions ####
#################################

npas_max = 20 # Nombre d'iterations algo K-means
maxiter = 50 # Nombre d'iterations algo EM
epsilon = 10^{-3} # Seuil de précision algo EM

####################
#### Question 8 ####
####################
# Execution de l'algorithme K-means
Kmeans(data, K, npas_max, T)

# Distribution du nombre de pas pour K-means sur n essais
npas_Kmeans(data, K, npas_max, n=20)

####################
#### Question 9 ####
####################
Algo_EM(data, K, npas_max, maxiter, epsilon, afficher = T, afficher_cercle = T)

#####################
#### Question 10 ####
#####################
Algo_EM(data, K, npas_max, maxiter, epsilon, afficher_ellipse = T)

########################
#### Question bonus ####
########################
Algo_EM(data, K, npas_max, maxiter, epsilon, afficher = T, afficher_ellipse = T, afficher_cercle = T)

