#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"




double potential_energy(Slice *psl) { 

    Pts *psi,*psj;
    vector dr[sys.npart][sys.npart],rij,p[sys.npart][sys.nsites],rnorm,zvec,zcen[sys.npart],pi,pj,ni,nj,n[2];
    double r2,r6,r6inv,r12inv,r24inv,r48inv,r,cosijtheta,cositheta,cosjtheta,itheta2,jtheta2;
    double potential_energy,potential_energy_attP=0,potential_energy_attC=0,potential_energy_rep=0;
    double potential_energy_attP_d,potential_energy_attC_d,potential_energy_rep_d;
    double phi_gauss,phi_i,phi_j,ensite,bondsmat[sys.npart][sys.npart];
    double enrep[sys.npart][sys.npart],eniso[sys.npart][sys.npart],enpatch[sys.npart][sys.npart],dist[sys.npart][sys.npart];
    double enrepmaxclust=0, enisomaxclust=0, enpatchmaxclust=0;
    double enrepmaxclustop=0, enisomaxclustop=0, enpatchmaxclustop=0;
    int clustsize, nbonds[sys.npart],trimer[3],posangle,negangle,npartforplane;
    int ipart,jpart,kpart,isite,jsite,iclust;
    tensor rotmati,rotmatj;

    //reset_center(psl);
    psl->maxclustersize=0;
    psl->enmaxclust=0;
    psl->nbonds=0;
    psl->bonddidnotbreak=0;
    psl->misalignments=0;

    zvec.x=0.;
    zvec.y=0.;
    zvec.z=1.;

    //printf("in poten\n");
  
    for( ipart=0; ipart<sys.npart; ipart++ )  {
        nbonds[ipart]=0;
        psi = &psl->pts[ipart];
        rotmati = getrotmatrix(psi->q);
        matrix_x_vector(rotmati,zvec,zcen[ipart]);
        for( isite=0; isite<sys.nsites; isite++) {
            matrix_x_vector(rotmati,sys.site[isite],p[ipart][isite]);
        }
        for(jpart=0; jpart<sys.npart; jpart++) {
            enrep[ipart][jpart]=0;
            eniso[ipart][jpart]=0;
            enpatch[ipart][jpart]=0;
            dr[ipart][jpart]=nulvec;
            dist[ipart][jpart]=10000;
            bondsmat[ipart][jpart]=0;
        }
    }

    //calculate Lennard-Jones interaction between all particles
    for( ipart=1; ipart<sys.npart; ipart++ )  {
        psi = &psl->pts[ipart];
        for( jpart=0; jpart<ipart; jpart++ ) {
            psj = &psl->pts[jpart];
            //calculate isotropic potential
            vector_minus(psi->r,psj->r,rij);
            pbc(rij,sys.boxl);
            dr[jpart][ipart]=rij;
            scalar_times(rij,-1.,dr[ipart][jpart]);
            r2 = vector_inp(rij,rij);
            dist[ipart][jpart]=r2;
            dist[jpart][ipart]=r2;
            if(r2<sys.rcutoffsq) {
                r6 = r2*r2*r2;
                r6inv=1.0/r6;
                r12inv = r6inv*r6inv;
                r24inv = r12inv*r12inv;
                if(r2<sys.sigmaLJsq) {
                    potential_energy_rep_d = r24inv-r12inv+0.25; 
                    if(potential_energy_rep>1000) {
                        printf("Warning: ipart %d jpart %d energy %lf distance %lf\n",
                                ipart,jpart,potential_energy_rep,r2);
                    }
                    enrep[ipart][jpart]=potential_energy_rep_d;
                    enrep[jpart][ipart]=potential_energy_rep_d;
                    potential_energy_rep += potential_energy_rep_d; 
                }
                potential_energy_attC_d = r24inv-r12inv; 
                potential_energy_attC += potential_energy_attC_d;
                eniso[ipart][jpart] = potential_energy_attC_d;
                eniso[jpart][ipart] = potential_energy_attC_d;
                r=sqrt(r2);
                scalar_divide(rij,r,rnorm);
                //calculate angular part of potential
                for( isite=0; isite<sys.nsites; isite++ ) {
                    cositheta=-vector_inp(rnorm,p[ipart][isite]);
                    //itheta2=acos(cositheta);
                    //itheta2*=itheta2;
                    itheta2 = sys.oneover72*(1.-cositheta)*(13.-cositheta)*(13.-cositheta);
                    //if(cositheta<sys.cosdelta) {
                    //    continue;
                    //}
                    //phi_i=0.5*(1.0-cos(PI*(cositheta-sys.cosdelta)*sys.oneover_cosdelta));
                    for( jsite=0; jsite<sys.nsites; jsite++ ) {
                        cosjtheta=vector_inp(rnorm,p[jpart][jsite]);
                        //jtheta2=acos(cosjtheta);
                        //jtheta2*=jtheta2;
                        jtheta2 = sys.oneover72*(1.-cosjtheta)*(13.-cosjtheta)*(13.-cosjtheta);
                        //if(cosjtheta<sys.cosdelta) {
                        //    continue;
                        //}
                        //phi_j=0.5*(1.0-cos(PI*(cosjtheta-sys.cosdelta)*sys.oneover_cosdelta));
                        phi_gauss = exp(-(itheta2+jtheta2)*sys.deltasq);
                        //ensite = phi_i*phi_j*(r24inv - r12inv);
                        ensite = phi_gauss*(r24inv - r12inv);
                        enpatch[ipart][jpart] += ensite;
                        enpatch[jpart][ipart] += ensite;
                        //printf("ipart %d jpart %d ensite %lf enpatch %lf\n",ipart,jpart,ensite,enpatch[ipart][jpart]);
                        //printf("ipart %d jpart %d isite %d jsite %d ensite %lf itheta2 %lf jtheta2 %lf phi_gauss %lf\n",ipart,jpart,isite,jsite,ensite,itheta2,jtheta2,phi_gauss);
                        //if(ensite < path.enbondclust) {
                        if(ensite < path.enbond) {
                            bondsmat[ipart][jpart]=1;
                            bondsmat[jpart][ipart]=1;
                            nbonds[ipart]++;
                            nbonds[jpart]++;
                        }
                        //if(ensite < path.enbond) {
                        //    nbonds[ipart]++;
                        //    nbonds[jpart]++;
                        //}
                        else if(ensite<path.enbondbreak) {
                            psl->bonddidnotbreak=1;
                        }
                        potential_energy_attP += ensite;
                    }
                }
            }
        }
    }

    int ov,label[sys.npart],notallgreysgone;
    //,clusterid[sys.npart]
    int nwhite, ngrey, nblack, ncluster, maxsize, WHITE=0, GREY=1, BLACK=2;
    double npart;
    vector rcom;

    //printf("performing cluster analysis\n");
    //do cluster analysis
    for(ipart=0; ipart<sys.npart; ipart++) {
        label[ipart]=WHITE;
        psl->clusterid[ipart]=0;
    }
    npart=0.;
    ncluster=0;
    nwhite = sys.npart;
    ngrey=0;
    nblack=0;
    for(ipart=0; ipart<sys.npart; ipart++) {
        if(label[ipart]==WHITE) {
            npart = 0.;
            ncluster++;
            nwhite--;
            ngrey++;
            label[ipart]=GREY;
            rcom=nulvec;
            do {
                notallgreysgone=0;
                for(jpart=0; jpart<sys.npart; jpart++) {
                    if(label[jpart]==GREY) {
                        notallgreysgone=1;
                        ngrey--;
                        nblack++;
                        label[jpart]=BLACK;
                        psl->clusterid[jpart]=ncluster;
                        //vector_add(psl->pts[jpart].r,rcom,rcom);
                        npart+=1.;
                        for(kpart=0; kpart<sys.npart; kpart++) {
                            if(kpart!=jpart) {
                                if(label[kpart]==WHITE) {   
                                    ov = bondsmat[jpart][kpart];
                                    if(ov==1) {
                                        label[kpart]=GREY;
                                        nwhite--;
                                        ngrey++;
                                    }
                                }
                            }
                        }
                    }
                }
                if(ngrey==0) notallgreysgone=0;
            } while(notallgreysgone);
            if(nwhite==0) {
                ipart=sys.npart;
            }
        }
    }


    //if(ncluster==1) {
    //    psl->mindist=0.;
    //}
    //else {
    //    psl->mindist=100000000000.;
    //}
    psl->mindist=100000000000.;
    for(ipart=0; ipart<sys.npart; ipart++) {
        for(jpart=0; jpart<sys.npart; jpart++) {
            if(ipart!=jpart) {
                if(psl->clusterid[ipart]!=psl->clusterid[jpart]) {
                    if(dist[ipart][jpart]<psl->mindist) {
                        psl->mindist=dist[ipart][jpart];
                        //printf("dist %d %d %lf\n",ipart,jpart,dist[ipart][jpart]);
                    }
                }
            }
        }
    }


    psl->maxclustersize=0;
    psl->indexmaxcluster=0;
    if(ncluster==1) {
        psl->maxclustersize=sys.npart;
        psl->indexmaxcluster=1;
    }
    else if (ncluster>1) {
        for(iclust=1; iclust<(ncluster+1); iclust++) {
            clustsize=0;
            for(jpart=0; jpart<sys.npart; jpart++) {
                if(psl->clusterid[jpart]==iclust) {
                    clustsize++;
                }
            }
            if(clustsize>psl->maxclustersize) {
                psl->maxclustersize=clustsize;
                psl->indexmaxcluster=iclust;
            }
        }
    }

    enrepmaxclust=0;
    enisomaxclust=0;
    enpatchmaxclust=0;
    enrepmaxclustop=0;
    enisomaxclustop=0;
    enpatchmaxclustop=0;
    for(ipart=0; ipart<sys.npart; ipart++) {
        if (psl->clusterid[ipart]==psl->indexmaxcluster) {
            for(jpart=ipart+1; jpart<sys.npart; jpart++) {
                if(psl->clusterid[jpart]==psl->indexmaxcluster) {
                    //printf("ipart %d jpart %d en %lf \n",ipart,jpart,enpatch[ipart][jpart]);
                    if(enpatch[ipart][jpart]<path.enbond) {
                        psl->nbonds++;
                    }

                    enrepmaxclust+=enrep[ipart][jpart];
                    enisomaxclust+=eniso[ipart][jpart];
                    enpatchmaxclust+=enpatch[ipart][jpart];
                }
            }
            //check for misalignments if particles have more than 1 bond
            if(nbonds[ipart]==2) {
                for(jpart=0; jpart<3; jpart++) {
                    trimer[0]=ipart;
                }
                trimer[0]=ipart;
                npartforplane=0;
                posangle=0;
                negangle=0;
	            for(jpart=0; jpart<sys.npart; jpart++) {
                    if(jpart!=ipart) {
                        if((bondsmat[ipart][jpart]==1) && (nbonds[jpart]>1)) {
                            trimer[npartforplane+1]=jpart;
                            if(psl->clusterid[jpart]==psl->indexmaxcluster) {
                                rij=dr[ipart][jpart];
                                r2=vector_inp(rij,rij);
                                scalar_divide(rij,sqrt(r2),rij);
                                n[npartforplane]=rij;
                                npartforplane++;
                                if(npartforplane>2) {
                                    printf("WARNING: something went wrong with plane calculation\n");
                                }
                            }
                        }
                    }
                }
                if(npartforplane==2) {
                    //try to not recount misalignments...
                    if( (nbonds[trimer[1]]==2) && (trimer[1]<trimer[0]) ) {
                        if( (nbonds[trimer[2]]==2) && (trimer[2]<trimer[0]) ) {
                            continue;
                        }
                    }
                    vector_cross(n[0],n[1],ni);
	                for(jpart=0; jpart<3; jpart++) {
                        cosijtheta=vector_inp(zcen[trimer[jpart]],ni);
                        if(cosijtheta>0) {
                            posangle++;
                        }
                        else if(cosijtheta<=0) {
                            negangle++;
                        }
	                }
                    if((posangle==3)||(negangle==3)) {
                        psl->misalignments+=0;
                        //continue;
                    }
                    else {
                        psl->misalignments++;
                    }
                }
            }
        }
        
        if (path.current_maxclustid[ipart]==1) {
            for(jpart=ipart+1; jpart<sys.npart; jpart++) {
                if (path.current_maxclustid[jpart]==1) {
                    enrepmaxclustop+=enrep[ipart][jpart];
                    enisomaxclustop+=eniso[ipart][jpart];
                    enpatchmaxclustop+=enpatch[ipart][jpart];
                }
            }
        }
    }
    psl->enmaxclust = 4.0*enrepmaxclust + 4.0*sys.epsilonP*enpatchmaxclust + 4.0*sys.epsilonC*enisomaxclust;
    psl->enmaxclustop = 4.0*enrepmaxclustop + 4.0*sys.epsilonP*enpatchmaxclustop + 4.0*sys.epsilonC*enisomaxclustop;

    potential_energy_rep*=4.0;
    potential_energy_attP*=4.0*sys.epsilonP;
    potential_energy_attC*=4.0*sys.epsilonC;
    potential_energy = potential_energy_rep + potential_energy_attC + potential_energy_attP;

    //printf("out poten\n");
    return potential_energy;
}



void calculate_forces(Slice *psl) {

    Pts *psi,*psj;
    vector rij,p[sys.npart][sys.nsites],ui,uj,rnorm;
    vector rcrosspi,rcrosspj,piperpr,pjperpr,fangi,fangj;
    double r2,r6,r6inv,r12inv,r24inv,r48inv,r2inv,rinv,potential_energy=0,potential_energy_rep=0,r,cositheta,cosjtheta;
    double phi_i,phi_j,fmag,fmagP,fmagrinv,UmagP;
    int ipart,jpart,isite,jsite;
    tensor rotmati,rotmatj;
    
    //initialize all forces to zero

    //printf("in calf force\n");
    for( ipart=0; ipart<sys.npart; ipart++ )  {
        psi = &psl->pts[ipart];
        psi->f = nulvec;
        psi->t = nulvec;
        rotmati = getrotmatrix(psi->q);
        for( isite=0; isite<sys.nsites; isite++) {
            matrix_x_vector(rotmati,sys.site[isite],p[ipart][isite]);
        }
    }

    //calculate forces based on Lennard-Jones interaction for all particles
    //add torques
    for( ipart=1; ipart<sys.npart; ipart++ ) {
        psi =&psl->pts[ipart];
        for( jpart=0; jpart<ipart; jpart++ ) {
            psj = &psl->pts[jpart];
            vector_minus(psi->r,psj->r,rij);
            pbc(rij,sys.boxl);
            r2 = vector_inp(rij,rij);
            if(r2<sys.rcutoffsq) {
                r2inv=1.0/r2;
                r6inv=r2inv*r2inv*r2inv;
                r12inv = r6inv*r6inv;
                r24inv = r12inv*r12inv;
                //r48inv = r24inv*r24inv;
                
                if(r2<sys.sigmaLJsq) {
                    fmag = r2inv*(96.*r24inv - 48.*r12inv);
                    //fmag = r2inv*(192.*r48inv - 96.*r24inv);
                    scalar_plustimes(rij,fmag,psi->f);
                    scalar_mintimes(rij,fmag,psj->f);
                }

                //fmag = sys.epsilonC*r2inv*(96.*r24inv - 48.*r12inv);
                ////fmag = sys.epsilonC*r2inv*(192.*r48inv - 96.*r24inv);
                //scalar_plustimes(rij,fmag,psi->f);
                //scalar_mintimes(rij,fmag,psj->f);

                UmagP = 4.0*sys.epsilonP*(r24inv - r12inv);
                //UmagP = 4.0*sys.epsilonP*(r48inv - r24inv);
                fmag = sys.epsilonP*r2inv*(96.*r24inv - 48.*r12inv);
                //fmag = sys.epsilonP*r2inv*(192.*r48inv - 96.*r24inv);

                r=sqrt(r2);
                //printf("distance particles %lf\n",r);
                rinv=1.0/r;
                scalar_times(rij,rinv,rnorm);
                //calculate angular part of force
                for( isite=0; isite<sys.nsites; isite++ ) {
                    cositheta=-vector_inp(rnorm,p[ipart][isite]);
                    if(cositheta<sys.cosdelta) {
                        continue;
                    }
                    //calculate phi for particle i
                    phi_i=0.5*(1.0-cos(PI*(cositheta-sys.cosdelta)*sys.oneover_cosdelta));
                    //printf("phi_i %d %lf\n",ipart, phi_i);
                    
                    for( jsite=0; jsite<sys.nsites; jsite++ ) {
                        cosjtheta=vector_inp(rnorm,p[jpart][jsite]);
                        if(cosjtheta<sys.cosdelta) {
                            continue;
                        }
                        //calculate phi for particle j
                        phi_j=0.5*(1.0-cos(PI*(cosjtheta-sys.cosdelta)*sys.oneover_cosdelta));
                        //printf("phi_j %d %lf\n",jpart, phi_j);

                        fmagP = fmag*phi_i*phi_j;
                        //printf("fmagP attraction %lf\n", fmagP);
                        scalar_plustimes(rij,fmagP,psi->f);
                        scalar_mintimes(rij,fmagP,psj->f);

                        vector_cross(rnorm,p[ipart][isite],rcrosspi);
                        vector_cross(rnorm,p[jpart][jsite],rcrosspj);

                        vector_cross(rnorm,rcrosspi,piperpr);
                        vector_cross(rnorm,rcrosspj,pjperpr);

                        //CALCULATE DEL U / DEL COSTHETA
                        fmag = UmagP*phi_j*HALFPI*sin(PI*(cositheta-sys.cosdelta)*sys.oneover_cosdelta)*sys.oneover_cosdelta; 
                        //printf("fmag due to patch i %lf\n", fmag);
                        fmagrinv=fmag*rinv;
                        scalar_mintimes(piperpr,fmagrinv,psi->f);
                        scalar_plustimes(piperpr,fmagrinv,psj->f);
                        //scalar_plustimes(piperpr,fmagrinv,psi->f);
                        //scalar_mintimes(piperpr,fmagrinv,psj->f);

                        scalar_mintimes(rcrosspi,fmag,psi->t);

                        fmag = UmagP*phi_i*HALFPI*sin(PI*(cosjtheta-sys.cosdelta)*sys.oneover_cosdelta)*sys.oneover_cosdelta; 
                        //printf("fmag due to patch j %lf\n", fmag);
                        fmagrinv=fmag*rinv;
                        scalar_plustimes(pjperpr,fmagrinv,psi->f);
                        scalar_mintimes(pjperpr,fmagrinv,psj->f);

                        scalar_plustimes(rcrosspj,fmag,psj->t);
                        
                        //printf("torque %d %lf %lf %lf\n",ipart,psi->t.x,psi->t.y,psi->t.z);
                        //printf("torque %d %lf %lf %lf\n",jpart,psj->t.x,psj->t.y,psj->t.z);
                    }
                }
            }
        }
    }

    //exit(1);

    //printf("out calf force\n");
    return ;
}



void propagate_bd(Slice *psl)  {

    int ipart,istep;
    Pts *psi;
    vector theta,f,t,u;
    quaternion qu1,qu2,qu3,qprime;
    tensor rotmat;
    quattensor Bmat;
    double dt,lambdaq;

    for(istep=0; istep<langevin.ninter; istep++) {
        calculate_forces(psl);

        for( ipart=0; ipart<sys.npart; ipart++) {
            psi = &psl->pts[ipart];

            //translational part
            //muT*F*dt*Beta
            scalar_times(psi->f,langevin.dtBeta,f);
            scalar_times(f,sys.mobilityT,f);
            vector_add(psi->r,f,psi->r);

            theta=RandomBrownianVector(langevin.dtD);
            scalar_times(theta,sys.sqrtmobilityT,theta);
            vector_add(psi->r,theta,psi->r);

            pbc(psi->r,sys.boxl);

            //rotational part
            //get matrices for q(t)
            rotmat = getrotmatrix(psi->q); 
            Bmat = getquatmatrix(psi->q);

            //Baalpha * (muR) * rotmatA * torque * delt = qu1
            scalar_times(psi->t,langevin.dtBeta,t);
            //do not forget, multiply with the inverse rotation matrix to convert to body-fixed torque
            matrixT_x_vector(rotmat,t,u);
            scalar_times(u,sys.mobilityR,u);
            quatmatrix_x_vec(Bmat,u,qu1);

            //Baalpha * (muR) * theta = qu2
            theta=RandomBrownianVector(langevin.dtD);
            scalar_times(theta,sys.sqrtmobilityR,theta);
            quatmatrix_x_vec(Bmat,theta,qu2);

            //qprime = qu1+qu2+q(t) for the first time
            quat_add(qu1,qu2,qprime);

            quat_add(qprime,psi->q,qprime);

            //find lambdaq, I guess it is th smallest...check which one keeps q(t+delt)^2=1
            //woohoo it works for the smaller lambdaq, maybe also simply for the bigger...
            lambdaq=langrange_multiplier_quat(qprime, psi->q);
            if(lambdaq>1e20) {
                //langrange multiplier did not work, simply renormalize quaternion
                scdivide_quat(qprime,sqrt(quat_inp(qprime,qprime)),psi->q);
            }
            else {
                //lambdaq*q(t) = qu3
                sctimes_quat(psi->q,lambdaq,qu3);
                //q(t+delt) = qprime+qu3
                quat_add(qprime,qu3,psi->q);
            }

        }
    }

    return;
}



void calculate_forces_nonoise(Slice *psl) {

    Pts *psi,*psj;
    vector rij,p[sys.npart][sys.nsites],ui,uj,rnorm;
    vector rcrosspi,rcrosspj,piperpr,pjperpr,fangi,fangj;
    double r2,r6,r6inv,r12inv,r24inv,r48inv,r2inv,rinv,potential_energy=0,potential_energy_rep=0,r,cositheta,cosjtheta;
    double phi_i,phi_j,fmag,fmagP,fmagrinv,UmagP;
    int ipart,jpart,isite,jsite;
    tensor rotmati,rotmatj;
    
    //initialize all forces to zero

    //printf("in calf force\n");
    for( ipart=0; ipart<sys.npart; ipart++ )  {
        if(psl->clusterid[ipart]==psl->indexmaxcluster) {
            psi = &psl->pts[ipart];
            psi->f = nulvec;
            psi->t = nulvec;
            rotmati = getrotmatrix(psi->q);
            for( isite=0; isite<sys.nsites; isite++) {
                matrix_x_vector(rotmati,sys.site[isite],p[ipart][isite]);
            }
        }
    }

    //calculate forces based on Lennard-Jones interaction for all particles
    //add torques
    for( ipart=1; ipart<sys.npart; ipart++ ) {
        if(psl->clusterid[ipart]==psl->indexmaxcluster) {
            psi =&psl->pts[ipart];
            for( jpart=0; jpart<ipart; jpart++ ) {
                if(psl->clusterid[jpart]==psl->indexmaxcluster) {
                    psj = &psl->pts[jpart];
                    vector_minus(psi->r,psj->r,rij);
                    pbc(rij,sys.boxl);
                    r2 = vector_inp(rij,rij);
                    if(r2<sys.rcutoffsq) {
                        r2inv=1.0/r2;
                        r6inv=r2inv*r2inv*r2inv;
                        r12inv = r6inv*r6inv;
                        r24inv = r12inv*r12inv;
                        //r48inv = r24inv*r24inv;
                        
                        if(r2<sys.sigmaLJsq) {
                            fmag = r2inv*(96.*r24inv - 48.*r12inv);
                            //fmag = r2inv*(192.*r48inv - 96.*r24inv);
                            scalar_plustimes(rij,fmag,psi->f);
                            scalar_mintimes(rij,fmag,psj->f);
                        }

                        fmag = sys.epsilonC*r2inv*(96.*r24inv - 48.*r12inv);
                        //fmag = sys.epsilonC*r2inv*(192.*r48inv - 96.*r24inv);
                        scalar_plustimes(rij,fmag,psi->f);
                        scalar_mintimes(rij,fmag,psj->f);

                        UmagP = 4.0*sys.epsilonP*(r24inv - r12inv);
                        //UmagP = 4.0*sys.epsilonP*(r48inv - r24inv);
                        fmag = sys.epsilonP*r2inv*(96.*r24inv - 48.*r12inv);
                        //fmag = sys.epsilonP*r2inv*(192.*r48inv - 96.*r24inv);

                        r=sqrt(r2);
                        //printf("distance particles %lf\n",r);
                        rinv=1.0/r;
                        scalar_times(rij,rinv,rnorm);
                        //calculate angular part of force
                        for( isite=0; isite<sys.nsites; isite++ ) {
                            cositheta=-vector_inp(rnorm,p[ipart][isite]);
                            if(cositheta<sys.cosdelta) {
                                continue;
                            }
                            //calculate phi for particle i
                            phi_i=0.5*(1.0-cos(PI*(cositheta-sys.cosdelta)*sys.oneover_cosdelta));
                            //printf("phi_i %d %lf\n",ipart, phi_i);
                            
                            for( jsite=0; jsite<sys.nsites; jsite++ ) {
                                cosjtheta=vector_inp(rnorm,p[jpart][jsite]);
                                if(cosjtheta<sys.cosdelta) {
                                    continue;
                                }
                                //calculate phi for particle j
                                phi_j=0.5*(1.0-cos(PI*(cosjtheta-sys.cosdelta)*sys.oneover_cosdelta));
                                //printf("phi_j %d %lf\n",jpart, phi_j);

                                fmagP = fmag*phi_i*phi_j;
                                //printf("fmagP attraction %lf\n", fmagP);
                                scalar_plustimes(rij,fmagP,psi->f);
                                scalar_mintimes(rij,fmagP,psj->f);

                                vector_cross(rnorm,p[ipart][isite],rcrosspi);
                                vector_cross(rnorm,p[jpart][jsite],rcrosspj);

                                vector_cross(rnorm,rcrosspi,piperpr);
                                vector_cross(rnorm,rcrosspj,pjperpr);

                                //CALCULATE DEL U / DEL COSTHETA
                                fmag = UmagP*phi_j*HALFPI*sin(PI*(cositheta-sys.cosdelta)*sys.oneover_cosdelta)*sys.oneover_cosdelta; 
                                //printf("fmag due to patch i %lf\n", fmag);
                                fmagrinv=fmag*rinv;
                                scalar_mintimes(piperpr,fmagrinv,psi->f);
                                scalar_plustimes(piperpr,fmagrinv,psj->f);
                                //scalar_plustimes(piperpr,fmagrinv,psi->f);
                                //scalar_mintimes(piperpr,fmagrinv,psj->f);

                                scalar_mintimes(rcrosspi,fmag,psi->t);

                                fmag = UmagP*phi_i*HALFPI*sin(PI*(cosjtheta-sys.cosdelta)*sys.oneover_cosdelta)*sys.oneover_cosdelta; 
                                //printf("fmag due to patch j %lf\n", fmag);
                                fmagrinv=fmag*rinv;
                                scalar_plustimes(pjperpr,fmagrinv,psi->f);
                                scalar_mintimes(pjperpr,fmagrinv,psj->f);

                                scalar_plustimes(rcrosspj,fmag,psj->t);
                                
                                //printf("torque %d %lf %lf %lf\n",ipart,psi->t.x,psi->t.y,psi->t.z);
                                //printf("torque %d %lf %lf %lf\n",jpart,psj->t.x,psj->t.y,psj->t.z);
                            }
                        }
                    }
                }
            }
        }
    }

    //exit(1);

    //printf("out calf force\n");
    return ;
}




void propagate_bd_nonoise(Slice *psl)  {

    int ipart,istep;
    Pts *psi;
    vector theta,f,t,u;
    quaternion qu1,qu2,qu3,qprime;
    tensor rotmat;
    quattensor Bmat;
    double dt,lambdaq;

    for(istep=0; istep<langevin.ninter; istep++) {
        calculate_forces_nonoise(psl);

        for( ipart=0; ipart<sys.npart; ipart++) {
            if(psl->clusterid[ipart]==psl->indexmaxcluster) {
                psi = &psl->pts[ipart];

                //translational part
                //muT*F*dt*Beta
                scalar_times(psi->f,langevin.dtBeta,f);
                scalar_times(f,sys.mobilityT,f);
                vector_add(psi->r,f,psi->r);


                pbc(psi->r,sys.boxl);

                //rotational part
                //get matrices for q(t)
                rotmat = getrotmatrix(psi->q); 
                Bmat = getquatmatrix(psi->q);

                //Baalpha * (muR) * rotmatA * torque * delt = qu1
                scalar_times(psi->t,langevin.dtBeta,t);
                //do not forget, multiply with the inverse rotation matrix to convert to body-fixed torque
                matrixT_x_vector(rotmat,t,u);
                scalar_times(u,sys.mobilityR,u);
                quatmatrix_x_vec(Bmat,u,qu1);


                //qprime = qu1+qu2+q(t) for the first time
                qprime=qu1;

                quat_add(qprime,psi->q,qprime);

                //find lambdaq, I guess it is th smallest...check which one keeps q(t+delt)^2=1
                //woohoo it works for the smaller lambdaq, maybe also simply for the bigger...
                lambdaq=langrange_multiplier_quat(qprime, psi->q);
                if(lambdaq>1e20) {
                    //langrange multiplier did not work, simply renormalize quaternion
                    scdivide_quat(qprime,sqrt(quat_inp(qprime,qprime)),psi->q);
                }
                else {
                    //lambdaq*q(t) = qu3
                    sctimes_quat(psi->q,lambdaq,qu3);
                    //q(t+delt) = qprime+qu3
                    quat_add(qprime,qu3,psi->q);
                }
            }
        }
    }

    return;
}


void find_minimum(Slice *psl) {

    Slice oldslice;
    double den,oldcosdelta,oldoneovercosdelta;
    //vector dr;
    //double dr2;
    //int ipart;

    oldcosdelta=sys.cosdelta;
    oldoneovercosdelta=sys.oneover_cosdelta;
    sys.cosdelta=cos(40.*(PI/180.));
    sys.oneover_cosdelta = 1.0/(1.0-sys.cosdelta);
    do {
        oldslice = *psl;
        propagate_bd_nonoise(psl);
        psl->energy = potential_energy(psl);
        //den = fabs(psl->energy - oldslice.energy);
        den = fabs(psl->enmaxclust - oldslice.enmaxclust);
    } while(den>1e-5);
    //} while(dr2>1e-15);

    sys.cosdelta = oldcosdelta;
    sys.oneover_cosdelta = oldoneovercosdelta;
    return;
}


//void find_minimum_forinstate(Slice *psl) {
//
//    Slice oldslice;
//    double den,oldcosdelta,oldoneovercosdelta,oldtimestep;
//    int i;
//    //int count=0;
//
//    oldcosdelta=sys.cosdelta;
//    oldoneovercosdelta=sys.oneover_cosdelta;
//    sys.cosdelta=cos(30.*(PI/180.));
//    sys.oneover_cosdelta = 1.0/(1.0-sys.cosdelta);
//    oldtimestep = langevin.dtBeta;
//    langevin.dtBeta*=3.0;
//    do {
//        //count++;
//        for(i=0; i<50; i++) {
//            oldslice = *psl;
//            propagate_bd_nonoise(psl);
//        }
//        psl->energy = potential_energy(psl);
//        den = fabs(psl->energy - oldslice.energy);
//    } while(den>1e-1);
//    //} while(count<100);
//
//    langevin.dtBeta = oldtimestep;
//    sys.cosdelta = oldcosdelta;
//    sys.oneover_cosdelta = oldoneovercosdelta;
//    return;
//}




tensor getrotmatrix(quaternion q) {

    tensor mat;
    double q0,q1,q2,q3,q0sq,q1sq,q2sq,q3sq;

    q0=q.q0;
    q1=q.q1;
    q2=q.q2;
    q3=q.q3;

    q0sq=q0*q0;
    q1sq=q1*q1;
    q2sq=q2*q2;
    q3sq=q3*q3;

    mat.x.x = q0sq + q1sq - q2sq - q3sq ;  mat.x.y = 2.0*(q1*q2 - q0*q3)       ;  mat.x.z = 2.0*(q1*q3 + q0*q2)       ;
    mat.y.x = 2.0*(q1*q2 + q0*q3)       ;  mat.y.y = q0sq - q1sq + q2sq - q3sq ;  mat.y.z = 2.0*(q2*q3 - q0*q1)       ;
    mat.z.x = 2.0*(q1*q3 - q0*q2)       ;  mat.z.y = 2.0*(q2*q3 + q0*q1)       ;  mat.z.z = q0sq - q1sq - q2sq + q3sq ; 

    return mat;
}




quattensor getquatmatrix(quaternion q) {

    quattensor qmat;
    double q0,q1,q2,q3;

    q0=0.5*q.q0;
    q1=0.5*q.q1;
    q2=0.5*q.q2;
    q3=0.5*q.q3;

    qmat.w.x = -q1 ;  qmat.w.y = -q2 ;  qmat.w.z = -q3 ;
    qmat.x.x =  q0 ;  qmat.x.y = -q3 ;  qmat.x.z =  q2 ;
    qmat.y.x =  q3 ;  qmat.y.y =  q0 ;  qmat.y.z = -q1 ;
    qmat.z.x = -q2 ;  qmat.z.y =  q1 ;  qmat.z.z =  q0 ;

    return qmat;
}




double langrange_multiplier_quat(quaternion qprime, quaternion q) {

    double lambdaq,a,b,c,det,root1,root2;

    //solving for equation 15 Ilie Briels den Otter
    //lambdaq^2 + 2*lambdaq*qprime*q + qprime*qprime == 1

    a = 1.0;
    b = 2.0*quat_inp(q,qprime);
    c = quat_inp(qprime,qprime)-1.0;

    det = b*b - 4.0*a*c;
    if(det<0.0) {
        //printf("Warning: determinant lower than 0, can not find roots\nWhat to do?\n");
        //printf("b value %lf c value %lf\n",b,c);
        return BIGNUM;
    }
    else {
        //hmmm which root to choose?
        root1=(-b + sqrt(det))/(2.0*a);
        root2=(-b + sqrt(det))/(2.0*a);
    }

    //go for minimum for now...makes more sense
    //if this does not work try just renormalizing it...
    
    if(fabs(root1)<fabs(root2)) {
        lambdaq=root1;
    }
    else {
        lambdaq=root2;
    }

    //printf("lambdaq %lf\n",lambdaq);
    return lambdaq;
}


void reset_center(Slice *psl) {

    int ipart;
    vector dr,rcom;

    for( ipart=1; ipart<sys.npart; ipart++) {
        vector_minus(psl->pts[ipart].r,psl->pts[0].r,dr);
        pbc(dr,sys.boxl);
        vector_add(dr,psl->pts[0].r,psl->pts[ipart].r);
    }

    rcom=nulvec;
    for(ipart=0; ipart<sys.npart; ipart++) {
        vector_add(rcom,psl->pts[ipart].r,rcom);
    }

    scalar_divide(rcom,sys.npart,rcom);
    for( ipart=0; ipart<sys.npart; ipart++) {
        vector_minus(psl->pts[ipart].r,rcom,psl->pts[ipart].r);
    }

    return;
}



//int clusteranalysis(Slice *psl) {
//
//    int ipart,jpart,kpart,jsite,ksite, ov, ngrey=0, nwhite=0, misalignments, nblack=0, npart, nbonds, ncluster, maxsize, notallgreysgone;
//    int label[sys.npart];
//    double energy;
//    Pts *psj,*psk;
//    int WHITE=0, GREY=1, BLACK=2;
//    vector boxl=sys.boxl;
//
//    for(ipart=0; ipart<sys.npart; ipart++) {
//        label[ipart]=WHITE;
//    }
//
//    nwhite = sys.npart;
//    maxsize=0;
//    npart=0;
//    nbonds=0;
//    ncluster=0;
//    misalignments=0;
//    nwhite--;
//    ngrey=1;
//    for(ipart=0; ipart<sys.npart; ipart++) {
//        if(label[ipart]==WHITE) {
//            npart = 0;
//            ncluster++;
//            nwhite--;
//            ngrey++;
//            label[ipart]=GREY;
//            do {
//                notallgreysgone=0;
//                for(jpart=0; jpart<sys.npart; jpart++) {
//                    if(label[jpart]==GREY) {
//                        notallgreysgone=1;
//                        ngrey--;
//                        nblack++;
//                        label[jpart]=BLACK;
//                        npart++;
//                        if(npart>maxsize) {
//                            maxsize=npart;
//                        }
//                        for(kpart=0; kpart<sys.npart; kpart++) {
//                            if(kpart!=jpart) {
//                                if(label[kpart]==WHITE) {   
//                                    ov = check_bond(&psl->pts[jpart],&psl->pts[kpart],&misalignments);
//                                    if(ov==1) {
//                                        label[kpart]=GREY;
//                                        nwhite--;
//                                        ngrey++;
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//                if(ngrey==0) notallgreysgone=0;
//            } while(notallgreysgone);
//            if(nwhite==0) {
//                ipart=sys.npart;
//            }
//        }
//    }
//
//    if(misalignments<=1) {
//        //printf("misalignments %d\n",misalignments);
//        return 0;
//    }
//    if(ncluster>=3) {
//        return 0;
//    }
//
//
//    return 1;
//}



