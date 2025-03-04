#include "_fusee.h"

//----------------------------------------------------------------------
// Détermination des variables nécessaires aux calculs

double pourcent_eau(fusee r){
    return r.water_volume / r.empty_volume;
}

double temp2(fusee r, env e){
    double res;
    res = (double)1-pourcent_eau(r);
    res = pow(res,e.gamma-1);
    return res*e.T_ext;
}

double beta(env e){
    return 1.03+(0.021*e.gamma);
}

double press2(fusee r, env e){
    double res;
    res = (r.empty_volume - r.water_volume)/r.empty_volume;
    res = pow(res,e.gamma);
    return r.pressure*res;
}

double p_trans(env e){
    double res;
    res = (double)(e.gamma+1)/2;
    res = pow(res,e.gamma/(e.gamma-1));
    return e.P_ext*res; 
}

double r_m = _R / 0.029;

double _c2(fusee r, env e){
    double t2 = temp2(r,e);
    double res;
    res = (double)e.gamma*t2*r_m;
    return sqrt(res);
}

//----------------------------------------------------------------------
// Détermination des paramètres de fin de la phase de propulsion

double imp(fusee r, env e){ // Renvoie l'augmentation de la vitesse de la fusee après le "gaz blast"
    double p2 = press2(r,e);
    double b = beta(e);
    double pt = p_trans(e);
    double c2 = _c2(r,e);
    

    double terme1 = p2/(b*e.P_ext); terme1 = pow(terme1, (e.gamma-1)/(2*e.gamma)); terme1 = 1-terme1;
    terme1 = pt*terme1; terme1 = terme1 / (p2*(e.gamma-1));
    double terme2 = b*e.P_ext/p2; terme2 = pow(terme2, (e.gamma+1)/(2*e.gamma)); terme2 = 1-terme2;
    double terme3 = 8/(e.gamma+1); terme3 = sqrt(terme3);
    double terme4 = p2*r.empty_volume/c2;
    double res = (double) terme3*terme4;
    res = res * (terme1 + terme2);
    return res/r.empty_mass;
}

double _tau(fusee r, env e){
    double c2 = _c2(r,e);
    double a_star = r.nozzle_radius*_PI*r.nozzle_radius;
    double res = (double) r.empty_volume/(a_star*c2);
    res = res*(2/(e.gamma-1));
    double terme1 = (e.gamma+1)/2;
    terme1 = pow(terme1, ((e.gamma+1)/(2*(e.gamma-1))));
    res = res*terme1;
    return res;
}

double t_gas(fusee r,env e){    // Calcul du temps d'expulsion du gaz
    double tau = _tau(r,e);
    double p2 = press2(r,e);
    double b = beta(e);
    double term1 = p2/(b*e.P_ext);
    term1 = pow(term1, ((e.gamma-1)/(2*e.gamma)));
    double res = (double)term1-1;
    return tau*res;
}

// 1ère méthode pour expulsion de l'eau :

double v_e(fusee r, env e){ // Vitesse d'expulsion de l'eau dans le goulot
    double p2 = press2(r,e);
    double res = (double) 2*(p2 - e.P_ext); // On prend p2 plutôt que r.pressure pour une meilleur approximation qui prend plus en compte le volume d'eau
    res = res/_RHO_EAU;
    return sqrt(res);
}

void expulsion_eau_v1(double T[3], fusee r, env e){ // Met a jour T avec les paramètres de la fin de l'expulsion de l'eau de la bouteille
    double ve = v_e(r,e);
    double a_star = r.nozzle_radius*_PI*r.nozzle_radius;
    double a = r.radius*_PI*r.radius;

    double debit = ve*a_star;
    
    double t_f = (double) r.water_volume/debit;    // Calcul du temps après expulsion de l'eau
    T[0] = t_f;

    double m0 = r.empty_mass + (_RHO_EAU*r.water_volume);

    double terme1 = (double)_RHO_EAU*a_star*ve*t_f; terme1 = terme1 / m0; terme1 = 1-terme1;    // Calcul de la vitesse après expulsion de l'eau
    double vf = (double) log(terme1);
    vf = -ve*vf; vf = vf - e.g*t_f;
    T[1] = vf;

    double terme2 = log(terme1);  // Calcul de la position après expulsion de l'eau
    double terme3 = m0/(_RHO_EAU*a);
    double zf = (double)ve*t_f*(1-terme2);
    zf += terme3*terme2;
    zf -= 0.5*e.g*pow(t_f,2);
    T[2] = zf;

}

void methode1(double T[3], fusee r, env e){ // T[0] = t_i ; T[1] = v_i ; T[2] = z_i, T contient les paramètres à la fin de la phase de propulsion
    if (r.pressure != 0 && r.water_volume != 0)
    {
        expulsion_eau_v1(T, r, e);
        double v = T[1];
        T[1] += imp(r,e);
        T[0] += t_gas(r,e);
        double v_moy = (v+ T[1])/2;
        T[2] += v_moy*t_gas(r,e);
    }
}


// 2ème méthode pour expulsion de l'eau :


void extrait_coefficients(char* polynome,long  double coefficients[10]) {
    for (int i = 0; i < 10; i+=1) {
        coefficients[i] = 0.0;
    }
    
    long double a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0, i = 0, j = 0;

    sscanf(polynome, "%Lf*x**9 %Lf*x**8 %Lf*x**7 %Lf*x**6 %Lf*x**5 %Lf*x**4 %Lf*x**3 %Lf*x**2 %Lf*x %Lf", &a, &b, &c, &d, &e, &f, &g, &h, &i, &j);
    
    coefficients[0] = a;
    coefficients[1] = b;
    coefficients[2] = c;
    coefficients[3] = d;
    coefficients[4] = e;
    coefficients[5] = f;
    coefficients[6] = g;
    coefficients[7] = h;
    coefficients[8] = i;
    coefficients[9] = j;
}


void methode2_false(fusee r){   // NE MARCHE PAS (sûrement à cause du sscanf de la fonction précédente) : D'après tests, toutes les cases de coef sont à 0.0 
    FILE *file = fopen("resultats.txt", "r");
    fseek(file, 0, SEEK_END);
    long length = ftell(file);
    fseek(file, 0, SEEK_SET);

    char *content = malloc(length + 1);

    fread(content, 1, length, file);
    content[length] = '\0';
    printf("\n%s\n", content);

    fclose(file);
    
    long double coef[10];
    extrait_coefficients(content, coef);
    free(content);
    
    double a_star = r.nozzle_radius*_PI*r.nozzle_radius;

    long double terme1 = (long double) _RHO_EAU*a_star;
    terme1 = terme1/((2*r.empty_mass + _RHO_EAU*r.water_volume)/2);
    long double v = 0.0;

    for (int i = 0; i < 10; i+=1)
    {
        v += coef[i]*pow(0.04,i); 
    }
    printf("\n\n zvizbv %Lf\n\n", v);
    
}



double methode2(fusee r, env e){  // Compromis valable pour 7 bar et 0.5L d'eau avec ces coefficients
    double T[3];
    expulsion_eau_v1(T,r,e);

    long double coef[10];
    coef[0] = 0;
    coef[1] = 1201.97708;
    coef[2] = 0-8606.05281;
    coef[3] = 78982.47317;
    coef[4] = -745940.98508;
    coef[5] = 5721645.37635;
    coef[8] = 0-884287017.88062;
    coef[7] = 191281156.60597;
    coef[6] = 0-32751295.64697;
    coef[9] = 191281156.60597;

    long double v = 0.0;

    for (int i = 0; i < 10; i+=1)
    {
        v += coef[i]*pow(0.04,i);
    }
    
    return (double)v;
}

void moy_methode(double T[3], fusee r, env e){  //Test pas réellement convaincant
    double v2 = methode2(r,e) + imp(r,e);
    methode1(T,r,e);
    T[1] = (T[1]+v2)/2;
}

// On remarque un manque de précision de la méthode n°2, un polynôme de dégré 4 est sûrement trop peu.