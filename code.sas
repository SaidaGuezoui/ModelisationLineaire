
/* Autheurs : Saida Guezoui et Benoit Gilles */

%let path =  C:\Users\saida\OneDrive\Desktop\Master 1\Modele lineaire;


/************************* Première analyse des données *************************/

/* Importation et corrélation */ 

DATA ozone ;
INFILE "&path\Tp2\ozone.csv"  FIRSTOBS=2 DELIMITER=";";
INPUT  obs maxO3 T9  T12  T15  Ne9 Ne12 Ne15 Vx9  Vx12 Vx15  maxO3v vent $ pluie $;
drop vent pluie ;
RUN;
Proc contents data = ozone ; 
run;

/*Analyse des corrélations */ 

PROC CORR DATA= ozone(drop=obs) 
		  plots = matrix(histogram); 
RUN;

 
/* Régression multivariée de MaxO3*/ 

PROC REG DATA= ozone ;
model maxO3=T9  T12  T15 Ne9 Ne12 Ne15 Vx9 Vx12  Vx15  maxO3v/ 
			CLB
			alpha=.5
			covb;
id obs ; 
RUN;


/************************* Prévision *************************/
 
/* Prédiction d'ozone */ 
proc reg data = ozone ; 
model maxO3= T9  T12  T15 Ne9 Ne12 Ne15 Vx9 Vx12  Vx15  maxO3v ; /* expliquer max03 par les 10 var*/ 
id obs ; 
output out = pred_ozone
	residual = residu 
	predicted = val_pred 
	lclm = Ic_inf_95 /* IC*/
	uclm= Ic_sup_95
	lcl= pred_b_95
	ucl= pred_h_95 
	press=press;
RUN; 

proc print data=pred_ozone (drop= obs T9  T12  T15 Ne9 Ne12 Ne15 Vx9 Vx12  Vx15  maxO3v);
run;

/* Calcul des résidus */
data pred_ozone ; 
set pred_ozone ; 
residu2 = residu * residu ; 
run; 

proc print data= work.pred_ozone;
run;

/* Calcul de la moyenne des résidus^2 (erreur d'ajustement)*/  
proc means data= pred_ozone n mean sum ; 
var residu2; 
run; 


/* Calcul des press^2*/ 

data pred_ozone; 
set pred_ozone; 
press2= press * press  ; 
run; 


/* Calcul de la moyenne des press ^2 (erreur de prédiction) */ 

proc means data = pred_ozone n mean sum ; 
var press2 ; 
run; 


/* Selection de variables avec la méthode RSQUARE */ 

PROC REG DATA=ozone plots (only) =rsquare cp; 
model maxO3=T9  T12  T15 Ne9 Ne12 Ne15 Vx9 Vx12  Vx15  maxO3v/ 
	  Selection= rsquare best=1 ; /*le meilleur parmi les modèles pour tous les crières*/ 
id obs;  
RUN;  


/* Selection de variables avec la méthode FORWARD */ 

PROC REG DATA=ozone plots (only) = none; 
model maxO3=T9  T12  T15 Ne9 Ne12 Ne15 Vx9 Vx12  Vx15  maxO3v/ 
	  Selection= forward
      slentry = 0.05 
      clb 
      alpha = .05 
      covb;   
id obs; 
output out= forward_ozone 
	residual = residu 
	predicted= prediction_moy 
	lclm = moyenne_basse_95
	uclm = moyenne_haute_95 
    lcl = pred_basse_95
    ucl = pred_haute_95 
	h= leverage 
	press = press ; 
RUN; 
 

/* Regression sur le sous modèle choisi*/ 

/*Corrélation*/ 
PROC CORR DATA= ozone(drop=obs T9 T15 Ne12 Ne15 Vx12 Vx15) 
		  plots = matrix(histogram); 
RUN;

PROC REG DATA= ozone ;
model maxO3 = T12 Ne9 Vx9 maxO3v/ 
			CLB
			alpha=.5
			covb;
id obs ; 
RUN;

/*Calcul des press^2*/ 

data forward_ozone;
set forward_ozone; 
press2= press * press  ; 
run; 

/* Calcul de la moyenne des press^2 (erreur de prédiction après selection des var*/ 

proc means data = forward_ozone n mean sum ; 
var press2 ; 
run; 


/************************* Problème de colinéarité ******************************/

/*Matrice de corrélation */ 
proc corr data=ozone(drop=obs)
plots=matrix(histogram);
run;

/*Modèle de régression : Regression multivariée de MaxO3 */ 

proc reg data = ozone;
model MaxO3 = T9 T12 T15 Ne9 Ne12 Ne15 Vx9 Vx12 Vx15 maxO3v /
clb  /*IC paramètres*/
alpha=0.05  /*niveau 95%*/
covb;  /* cov estimateurs paramètres*/
id obs; /* identifiant de chaque obs*/
run;

 
/*Regression avec diagnostic de collinéarité par VIF*/ 

proc reg data = ozone plots=none ;
model MaxO3 = T9 T12 T15 Ne9 Ne12 Ne15 Vx9 Vx12 Vx15 maxO3v /
vif /* VIF */
tol ;/* valeur de tolérence pour les estimations */

id obs; /* identifiant de chaque obs*/
run;

 
/*Regression avec diagnostic de collinéarité par valeurs propres*/ 

proc reg data = ozone plots=none;
model MaxO3 = T9 T12 T15 Ne9 Ne12 Ne15 Vx9 Vx12 Vx15 maxO3v /
		collinoint; /* valeurs propres de X'X*/
id obs; /* identifiant de chaque obs*/
run;



/********************** Traitement de la multicolinéarité ********************/ 

 
/*Regression ridge, choix de pénalisation*/

proc reg data = ozone plots(only)=ridge
outest= b outvif ridge=0 to 0.1 by .005;
model MaxO3 = T9 T12 T15 Ne9 Ne12 Ne15 Vx9 Vx12 Vx15 maxO3v;
id obs; /* identifiant de chaque obs*/
run;


/*Regression ridge, pénalisation choisie*/

proc reg data = ozone outest=b2 plots=none;
model MaxO3 = T9 T12 T15 Ne9 Ne12 Ne15 Vx9 Vx12 Vx15 maxO3v /
ridge=0.015;
id obs; /* identifiant de chaque obs*/
run;


proc score data = ozone score = b2
out = ridgeprev predict type=ridge;
var T9 T12 T15 Ne9 Ne12 Ne15 Vx9 Vx12 Vx15 maxO3v;
id obs;
run;

data ridgeprev;
merge ozone ridgeprev;
by obs;
run;

data ridgeprev;
set ridgeprev;
residu = maxO3 - model1;
residu2 = residu * residu;
run;

proc means data = ridgeprev n mean;
var residu2;
run;

/* Régression sur composantes principales */

Proc princomp data = ozone out = composante; 
var T9 T12 T15 Ne9 Ne12 Ne15 Vx9 Vx12 Vx15 maxO3v;
run; 

proc reg data= composante;
model maxO3= prin1--prin10 / selection = rsquare cp best = 1; 
run; 

/*Regression du sous modèle */ 

PROC REG DATA= ozone ;
model maxO3= T12 Ne9 Vx9  maxO3v/ 
			CLB
			alpha=.5
			covb;
id obs ; 
RUN;
