#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "path.h"

//test for git


void setup_simulation() {

    int istate;
    printf("Setting up the system\n");

    InitializeRandomNumberGenerator(time(0l));

    read_input();
    init_model();

    printf("Allocating memory for the pathways\n");
    //for bmdcycle only initial positions need to be set
    if(sys.sim_type==0) {
        slice = (Slice *)calloc(1,sizeof(Slice));
        //setup_positions(&slice[0]);
        read_lambda();
        state_init();
        printf("Initial Energy system %lf\n",potential_energy(&slice[0]));
    }
    else if(sys.sim_type==1) {
        slice = (Slice *)calloc(MAXSLICES,sizeof(Slice));
        trial = (Slice *)calloc(MAXSLICES,sizeof(Slice));
        fwtrial = (Slice *)calloc(MAXSLICES,sizeof(Slice));
        bwtrial = (Slice *)calloc(MAXSLICES,sizeof(Slice));

        read_lambda();
        //target_input();
        state_init();
        replica_init();
        if(sys.start_type==1) {
            target_input();
            dos_input();
            traj_input();
        }
        stats_init();

        print_tis_setup();
    }

    printf("Setting up the system is done\n");

    return;
}


void init_model() {

    int isite,jsite;
    double angle,r2,l,psifactor,psi;
    vector u,zvec;
    tensor rotmat;
    quaternion qrot;

    sys.boxl.y = sys.boxl.z = sys.boxl.x;

    sys.oneover72 = 1./72.;

    sys.temp = 1.0/sys.beta;

    sys.rcutoffsq=sys.rcutoff*sys.rcutoff;
    printf("rcutoffsq %lf\n",sys.rcutoffsq);

    sys.sigmaLJ = pow(2.0,(1.0/12.0));
    //sys.sigmaLJ = pow(2.0,(1.0/24.0));
    printf("sigmaLJ   %lf\n",sys.sigmaLJ);
    sys.sigmaLJsq = sys.sigmaLJ*sys.sigmaLJ;
    printf("sigmaLJsq %lf\n",sys.sigmaLJsq);
    sys.delta=PI*sys.delta/180.0;
    sys.deltasq=1.0/(sys.delta*sys.delta);
    sys.cosdelta=cos(sys.delta);
    sys.oneover_cosdelta = 1.0/(1.0-sys.cosdelta);
    langevin.dtD=sqrt(2.0*langevin.timestep);
    langevin.dtBeta=langevin.timestep*sys.beta;

    sys.cosmisalignangle = cos(PI*sys.misalignangle/180.);
    printf("cos(misangle) %lf\n",sys.cosmisalignangle);

    sys.sqrtmobilityT = sqrt(sys.mobilityT);
    sys.sqrtmobilityR = sqrt(sys.mobilityR);

    if(sys.pol_type==1) {
        psifactor=sqrt(2.0/3.0);
        for(isite=0; isite<sys.nsites; isite++) {
            sys.site[isite].x=cos(isite*120.0*PI/180.0)*sqrt(1.0-psifactor*psifactor);
            sys.site[isite].y=sin(isite*120.0*PI/180.0)*sqrt(1.0-psifactor*psifactor);
            sys.site[isite].z=psifactor;
        }
    }

    else if(sys.pol_type==2) {
        psifactor = 0.5*sqrt(2.);
        for(isite=0; isite<sys.nsites; isite++) {
            sys.site[isite].x=cos(isite*90.0*PI/180.0)*sqrt(1.0-psifactor*psifactor);
            sys.site[isite].y=sin(isite*90.0*PI/180.0)*sqrt(1.0-psifactor*psifactor);
            sys.site[isite].z=psifactor;
        }
    }


    else if(sys.pol_type==3) {
        psifactor=0.3568220897730899;
        for(isite=0; isite<sys.nsites; isite++) {
            sys.site[isite].x=cos(isite*120.0*PI/180.0)*sqrt(1.0-psifactor*psifactor);
            sys.site[isite].y=sin(isite*120.0*PI/180.0)*sqrt(1.0-psifactor*psifactor);
            sys.site[isite].z=psifactor;
        }
    }
    else if(sys.pol_type==4) {
        psi=(1.+sqrt(5.))/2.0;
        psifactor = 1.0/((sqrt(psi*psi+1.0)));
        for(isite=0; isite<sys.nsites; isite++) {
            sys.site[isite].x=cos(isite*72.0*PI/180.0)*sqrt(1.0-psifactor*psifactor);
            sys.site[isite].y=sin(isite*72.0*PI/180.0)*sqrt(1.0-psifactor*psifactor);
            sys.site[isite].z=psifactor;
        }
    }
    else if(sys.pol_type==5) {
        sys.site[0].x = -1.;
        sys.site[0].y = 0.;
        sys.site[0].z = 0.;
        sys.site[1].x = 0.5;
        sys.site[1].y = -0.8090169943749475;
        sys.site[1].z = 0.3090169943749472;
        sys.site[2].x = 0.5;
        sys.site[2].y = 0.8090169943749475;
        sys.site[2].z = 0.3090169943749472;

        angle=5.820;
        angle+=0.00001;
        qrot.q0=cos(angle*(PI/180.));
        qrot.q1=0.0*sin(angle*(PI/180.));
        qrot.q2=1.0*sin(angle*(PI/180.));
        qrot.q3=0.0*sin(angle*(PI/180.));
        rotmat=getrotmatrix(qrot);
        //for(isite=0; isite<sys.nsites; isite++) {
        //    vprint(sys.site[isite]);
        //}
        for(isite=0; isite<sys.nsites; isite++) {
            matrix_x_vector(rotmat,sys.site[isite],u);
            sys.site[isite]=u;
            vprint(u);
        }
    }

    else if(sys.pol_type==6) {
        psifactor = sqrt(3.)/3.;
        for(isite=0; isite<sys.nsites; isite++) {
            sys.site[isite].x=cos(isite*120.0*PI/180.0)*sqrt(1.0-psifactor*psifactor);
            sys.site[isite].y=sin(isite*120.0*PI/180.0)*sqrt(1.0-psifactor*psifactor);
            sys.site[isite].z=psifactor;
            vprint(sys.site[isite]);
        }
    }


    for(isite=0; isite<sys.nsites; isite++) {
        printf("length patch vector %d: %lf\n",isite, vector_inp(sys.site[isite],sys.site[isite]));
    }
    for(isite=1; isite<sys.nsites; isite++) {
        for(jsite=0; jsite<isite; jsite++) {
            psi=vector_inp(sys.site[isite],sys.site[jsite]);
            printf("Angle between patch %d and %d: %lf\n",isite,jsite,(180./PI)*acos(psi));
        }
    }

    path.enbond=-0.65*0.25; //half of the max energy (note that the max energy: 4*eps)
    path.enbondclust=-0.65*0.25; //half of the max energy (note that the max energy: 4*eps)
    path.enbondbreak=-0.005*0.25; //half of the max energy (note that the max energy: 4*eps)
    //exit(1);

    return;
}




void setup_positions(Slice *psl) {

    double boxlength=sys.boxl.x, r2;
    int ipart, jpart, overlap;
    vector dr;
    Pts *psi, *psj;

    printf("Setting the positions for all the particles\n");
    //random configuration
    if(sys.npart==1) {
        printf("Only one particle: placing it at the centre\n");
        psi=&psl->pts[0];
        psi->r=nulvec;
        psi->q=RandomQuaternion();
    }

    else if(sys.npart==2) {
        printf("Only two particles: placing them towards each other\n");
        for( ipart=0; ipart<sys.npart; ipart++ ) {
            psi=&psl->pts[ipart];
            psi->r=nulvec;
            psi->r.z -= 0.5 - sys.sigmaLJ*ipart*sys.sigma;
            //psi->q=RandomQuaternion();
            if(ipart==0) {
                psi->q.q0 = 0.0;
                psi->q.q1 = 0.0;
                psi->q.q2 = 0.0;
                psi->q.q3 = 1.0;
            }
            if(ipart==1) {
                psi->q.q0 = 0.0;
                psi->q.q1 = 1.0;
                psi->q.q2 = 0.0;
                psi->q.q3 = 0.0;
            }
        }
    }


    else if(sys.start_type==1) {
        printf("Reading from conf.inp\n");
        conf_input(&slice[0]);
    }
    else if(sys.start_type==2) {
        if(sys.pol_type==1) {
            printf("Setting up tetrahedron\n");
            setup_tetrahedron(&slice[0]);
        }
        else if(sys.pol_type==2) {
            printf("Setting up tetrahedron\n");
            setup_octahedron(&slice[0]);
        }
        else if(sys.pol_type==3) {
            printf("Setting up dodecahedron\n");
            setup_dodecahedron(&slice[0]);
        }
        else if(sys.pol_type==4) {
            printf("Setting up icosahedron\n");
            setup_icosahedron(&slice[0]);
        }
        else if(sys.pol_type==5) {
            printf("Setting up truncated icosahedron\n");
            setup_truncatedicosahedron(&slice[0]);
        }
        else if(sys.pol_type==6) {
            printf("Setting up cube\n");
            setup_cube(&slice[0]);
        }
    }

    else if(sys.start_type==0) {
        printf("Randomly placing particles with random orientation\n");
        for( ipart=0; ipart<sys.npart; ipart++) {
            psi=&psl->pts[ipart];
            do {
                overlap=0;
                psi->r=RandomVector(boxlength);
                psi->q=RandomQuaternion();
                for( jpart=0; jpart<ipart; jpart++) {
                    psj=&psl->pts[jpart];
                    vector_minus(psj->r,psi->r,dr);
                    pbc(dr,sys.boxl);
                    r2=vector_inp(dr,dr);
                    if(r2<=sys.sigmaLJsq) {
                        overlap=1;
                    }
                }
            } while(overlap==1);
        }
    }

    printf("Setting the positions for all the particles is done\n");
    return;
}


#define MAXLINE 100
void read_input() {

    FILE *fp;
    char *pt,line[MAXLINE];

    if((fp = fopen("path.inp","r"))==NULL) {
        printf("ERROR: could not read input file\n");
        exit(1);
    }

    printf("\nReading input from path.inp\n");

    while(fgets(line,MAXLINE, fp) != NULL) {
        pt = strtok(line," ");
        if(strcmp(pt,"ncycle1")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.ncycle1);
        } else if( strcmp(pt,"ncycle2")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.ncycle2);
        } else if( strcmp(pt,"npart")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.npart);
        } else if( strcmp(pt,"graphics")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.graphics);
        } else if( strcmp(pt,"delta")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.delta);
        } else if( strcmp(pt,"diffenrep")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&path.diffenrep);
        } else if( strcmp(pt,"nsites")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.nsites);
        } else if( strcmp(pt,"beta")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.beta);
        } else if( strcmp(pt,"sigma")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.sigma);
        } else if( strcmp(pt,"rcutoff")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.rcutoff);
        } else if( strcmp(pt,"misalignangle")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.misalignangle);
        } else if( strcmp(pt,"mobilityT")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.mobilityT);
        } else if( strcmp(pt,"mobilityR")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.mobilityR);
        } else if( strcmp(pt,"epsilonC")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.epsilonC);
        } else if( strcmp(pt,"epsilonP")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.epsilonP);
        } else if( strcmp(pt,"boxl")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&sys.boxl.x);
        } else if( strcmp(pt,"timestep")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%lf",&langevin.timestep);
        } else if( strcmp(pt,"ninter")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&langevin.ninter);
        } else if( strcmp(pt,"pol_type")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.pol_type);
        } else if( strcmp(pt,"findnewstates")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.findnewstates);
        } else if( strcmp(pt,"sim_type")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.sim_type);
        } else if( strcmp(pt,"start_type")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&sys.start_type);
        } else if( strcmp(pt,"nshoot")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&path.nshoot);
        } else if( strcmp(pt,"nrepswap")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&path.nrepswap);
        } else if( strcmp(pt,"nstateswap")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&path.nswapstates);
        } else if( strcmp(pt,"nreverse")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&path.nreverse);
        } else if( strcmp(pt,"stateswapbias")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&path.stateswapbias);
        } else if( strcmp(pt,"fixedbias")==0) {
            pt = strtok(NULL," ");
            sscanf(pt,"%d",&path.fixedbias);
        } else if( strcmp(pt,"\n")==0) {
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"Simulation\n")==0) {
            printf("Reading simulation parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"Potential\n")==0) {
            printf("Reading potential parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"System\n")==0) {
            printf("Reading system parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"Path\n")==0) {
            printf("Reading path parameters\n");
            pt = strtok(NULL," ");
        } else if( strcmp(pt,"TIS\n")==0) {
            printf("Reading TIS parameters\n");
        } else {
            printf("Keyword unknown: %s\n",pt);
        }
    }

    fclose(fp);

    printf("Done reading path.inp\n");

    print_input();

    return;
}


void print_input(void) {

    printf("\nInput read from path.inp:\n");
    printf("Simulation\n");
    printf("sim_type        %d\n", sys.sim_type);
    printf("start_type      %d\n", sys.start_type);
    printf("pol_type        %d\n", sys.pol_type);
    printf("findnewstates   %d\n", sys.findnewstates);
    printf("graphics        %d\n", sys.graphics);
    printf("ncycle1         %d\n", sys.ncycle1);
    printf("ncycle2         %d\n", sys.ncycle2);

    printf("\n");
    printf("Potential\n");
    printf("sigma           %lf\n", sys.sigma);
    printf("epsilonC        %lf\n", sys.epsilonC);
    printf("epsilonP        %lf\n", sys.epsilonP);

    printf("\n");
    printf("System\n");
    printf("npart           %d\n", sys.npart);
    printf("nsites          %d\n", sys.nsites);
    printf("delta           %lf\n", sys.delta);
    printf("beta            %lf\n", sys.beta);
    printf("mobilityT       %lf\n", sys.mobilityT);
    printf("mobilityR       %lf\n", sys.mobilityR);
    printf("boxl            %lf\n", sys.boxl.x);
    printf("rcutoff         %lf\n", sys.rcutoff);
    printf("misalignangle   %lf\n", sys.misalignangle);

    printf("\n");
    printf("Path\n");
    printf("timestep        %lf\n", langevin.timestep);
    printf("ninter          %d\n", langevin.ninter);

    printf("\n");
    printf("TIS\n");
    printf("diffenrep       %lf\n",path.diffenrep);
    printf("nshoot          %d\n", path.nshoot);
    printf("nrepswap        %d\n", path.nrepswap);
    printf("nstateswap      %d\n", path.nswapstates);
    printf("nreverse        %d\n", path.nreverse);
    printf("stateswapbias   %d\n", path.stateswapbias);
    printf("fixedbias       %d\n", path.fixedbias);

    return;
}





void state_init() {

    int islice,ipart;

    //setup first slice as dodecahedron
    //  copy dodecahedron setup from tetrafull
    //  give each particle proper quaternion by rotating accordingly
    //
    //define the state volume boundaries
    //  unbound state has 0 energy, 0 bonds
    //  dodecahedron state has 30 bonds, therefore gs_en = -30*epsilonP energy
    //  trapped state has no specific state boundaries

    //set initial state at Icosahedron state


    path.enbond=-0.65*0.25; //half of the max energy (note that the max energy: 4*eps)
    path.enbondclust=-0.65*0.25; //half of the max energy (note that the max energy: 4*eps)
    path.enbondbreak=-0.005*0.25; //half of the max energy (note that the max energy: 4*eps)
    if(sys.pol_type==1) {
        setup_tetrahedron(&(slice[0]));
    }
    else if(sys.pol_type==2) {
        setup_octahedron(&(slice[0]));
    }
    else if(sys.pol_type==3) {
        setup_dodecahedron(&(slice[0]));
    }
    else if(sys.pol_type==4) {
        setup_icosahedron(&(slice[0]));
    }
    else if(sys.pol_type==5) {
        setup_truncatedicosahedron(&(slice[0]));
    }
    else if(sys.pol_type==6) {
        setup_cube(&(slice[0]));
    }

    slice[0].energy=potential_energy(&slice[0]);
    state[1].target=slice[0];
    state[1].target.mindist=sys.rcutoffsq;
    //keep the state as small as possible
    state[1].volume_op = state[1].lambda[0];
    printf("Definition of Icosahedron state:\n");
    printf("         Ground-state energy:     %lf\n",state[1].target.energy);
    printf("         Number of bonds:         %d\n",state[1].target.nbonds);
    printf("         State volume set at:     %lf\n",state[1].volume_op);
    printf("\n");

    state[0].target.energy=0;
    state[0].target.nbonds=0;
    state[0].target.mindist=sys.rcutoffsq;
    state[0].target.misalignments=0;
    state[0].target.maxclustersize=1;
    //keep the state as small as possible
    state[0].volume_op = state[0].lambda[0];
    printf("Definition of Unbound state:\n");
    printf("         Ground-state energy:     %lf\n",state[0].target.energy);
    printf("         Number of bonds:         %d\n",state[0].target.nbonds);
    printf("         State volume set at:     %lf\n",state[0].volume_op);
    printf("\n");


    path.nstates=2;
    path.initial_state=2;
    path.current_gsen = state[1].target.energy;

    //all particles belong to the maximum cluster initially
    for(ipart=0; ipart<sys.npart; ipart++) {
        path.current_maxclustid[ipart]=1;
    }

    slice[0] = state[1].target;
    create_all_rc(&(slice[0]));
    if(sys.sim_type==1) {
        for(islice=1; islice<MAXSLICES; islice++) {
            slice[islice]=slice[0];
        }
    }
    print_rc(&(slice[0]),path.initial_state);

    return;
}



void replica_init() {


    Replica *prep;
    int islice,istate,irep,jrep,pathlen,maxlength,len,type;

    //define interface volumes for every state
    //  interfaces for the unbound state are at very low energies
    //  interfaces for the Icosahedron state are starting from gsen(Ih) and up
    //  interfaces for trapped state are also starting from en min
    //
    
    
    printf("Replica definition per state:\n");
    for(istate=0; istate<MAXSTATES; istate++) {
        //printf("      State %d:\n",istate);
        //state[istate].nrep = state[1].nrep;
        path.nreplica = state[istate].nrep;
        //printf("            Number of replicas %d\n",path.nreplica);
        for( irep=0; irep<MAXREPLICA; irep++) {
            state[istate].srep[irep].index = irep;
            state[istate].srep[irep].swapindex = irep;
            //printf("            %d %lf\n",irep,state[istate].srep[irep].lambda);
        }
        for( irep=0; irep<path.nreplica; irep++) {
            if(istate<path.nstates) {
                state[istate].srep[irep].lambda = state[istate].lambda[irep];
                state[istate].lambda[irep] = state[istate].lambda[irep];
            }
        }
    }

    for(istate=0; istate<path.nstates; istate++) {
        printf("      State %d:\n",istate);
        path.nreplica = state[istate].nrep;
        printf("            Number of replicas %d\n",path.nreplica);
        for( irep=0; irep<path.nreplica; irep++) {
            printf("            %d %.10lf\n",irep,state[istate].srep[irep].lambda);
        }
    }


    path.nreplica = state[path.initial_state-1].nrep;
    for(irep=0; irep<path.nreplica; irep++) {
        replica[irep] = &state[path.initial_state-1].srep[irep];
    }


    path.scalefactor = 0.010;
    for(istate=0; istate<MAXSTATES; istate++) {
        state[istate].scalefactor = 0.010;
    }
    printf("Bootstrapping initial path for minus interface Icosahedron state\n");
    prep = replica[0];
    maxlength= MAXSLICES/2;
    len = trajectory_state_i(&slice[0],prep,maxlength,path.initial_state);

    printf("Length initial path from target to first interface: %d slices\n",len);

    for( islice=0 ;islice<=len; islice++ ) {
        slice[islice+len+1]=slice[islice];  
    }

    for( islice=0; islice<=len; islice++) {
        slice[len -islice] = slice[islice+len+1];
    }
  
    pathlen = 2*len+2;
    printf("Total length initial path minus interface: %d slices\n",pathlen);
    for( islice=0; islice<pathlen; islice++) {
        print_rc(&slice[islice],path.initial_state);
    }

    type =analyse_state_i(slice,prep,pathlen, path.initial_state);
    printf("Type initial path in minus interface: %d\n",type);

    replica[0]->type=type;
    replica[0]->pathlen= pathlen;
    path.nslices= pathlen;
    path.energy=slice[0].energy;
    path.current_replica=0;

    printf("Finished with replica initialization\n");

    return;
}



void read_lambda() {

    FILE *fplam;
    int irep;


    if((fplam = fopen("lambda.inp","r"))==NULL){
        printf("Warning: lambda.inp not found\n");
        exit(1);
    }
    fscanf(fplam,"%d",&state[1].nrep);
    for( irep=0; irep<state[1].nrep; irep++) {
        fscanf(fplam,"%lf",&state[1].lambda[irep]);
    }
    fclose(fplam);

    if((fplam = fopen("lambda_u.inp","r"))==NULL){
        printf("Warning: lambda_u.inp not found\n");
        exit(1);
    }
    fscanf(fplam,"%d",&state[0].nrep);
    for( irep=0; irep<state[0].nrep; irep++) {
        fscanf(fplam,"%lf",&state[0].lambda[irep]);
    }
    fclose(fplam);

    return;
}



void setup_tetrahedron(Slice *psl) {
    
    int ix,iy,iz,sign;
    int ipart,jpart,ipatch,isite;
    double phi,ru,scale;
    Pts *psi,*psj;
    double r2,l,sqrt2;
    vector dr;
    vector tetrasites[4][3],N,M,u0,u1;

    phi = (1.+sqrt(5.))/2.;

    printf("Setting all positions for the Tetrahedron and scaling bond lengths to sigma_LJ\n");
    sqrt2 = sqrt(2);

    ipart=0;
    for(ix=-1; ix<=1; ix+=2) {
        psl->pts[ipart].r.x = ix;
        psl->pts[ipart].r.z = -1/sqrt2;
        psl->pts[ipart].r.y = 0.;
        ipart++;
    }

    for(iy=-1; iy<=1; iy+=2) {
        psl->pts[ipart].r.y = iy;
        psl->pts[ipart].r.z = 1/sqrt2;
        psl->pts[ipart].r.x = 0.;
        ipart++;
    }


    ru = sqrt(3./8.)*sys.sigma;
    scale = sys.sigmaLJ*ru/sqrt(1.5);
    for(ipart=0; ipart<sys.npart; ipart++) {
        scalar_times(psl->pts[ipart].r,scale,psl->pts[ipart].r);
    }


    printf("Positions of all particles\n");
    for(ipart=0; ipart<sys.npart; ipart++) {
        printf("%lf %lf %lf\n",psl->pts[ipart].r.x,psl->pts[ipart].r.y,psl->pts[ipart].r.z);
    }



    printf("Obtaining correct quaternion for rotation\n");
    for(ipart=0; ipart<sys.npart; ipart++) {
        psi = &(psl->pts[ipart]);
        ipatch =0;
        for(jpart=0; jpart<sys.npart; jpart++) {
            if(jpart!=ipart) {
                psj = &(psl->pts[jpart]);
                vector_minus(psj->r,psi->r,dr);
                r2 = vector_inp(dr,dr);
                l = sqrt(r2);
                if(l<(sys.sigmaLJ+0.001)) {
                    scalar_divide(dr,l,dr);
                    tetrasites[ipart][ipatch] = dr;
                    //r2=vector_inp(tetrasites[ipart][ipatch],tetrasites[ipart][ipatch]);
                    //printf("length patch %d particle %d: %lf\n",ipatch,ipart,sqrt(r2));
                    ipatch++;
                }
            }
        }
        //dprint(ipatch);
    }

    //find proper quaternion for every particle
    //  first rotate vector (0,0,1) in proper direction of particle to centre of polyhedron
    //  subsequently rotate one of the patches along rotate (0,0,1) axis, other should immediately follow

    vector rotax;
    tensor rotmat;
    int foundquat,count;
    double angle;
    quaternion qrot,qrot2;
    vector u,zvec,zvecipart;

    zvec.x=0; zvec.y=0; zvec.z=1.0;

    for(ipart=0; ipart<sys.npart; ipart++) {
        //printf("Particle %d\n",ipart);
        psi = &(psl->pts[ipart]);
        scalar_divide(psi->r,-sqrt(vector_inp(psi->r,psi->r)),zvecipart);
        qrot = quatVecToVec(zvec,zvecipart);
        rotax = zvecipart;
        psi->q = qrot;
        angle = 0.001;
        qrot.q0 = cos(angle/2.0);
        qrot.q1 = rotax.x*sin(angle/2.0);
        qrot.q2 = rotax.y*sin(angle/2.0);
        qrot.q3 = rotax.z*sin(angle/2.0);
        count=0;
        do {
            foundquat=0;
            quat_times(qrot,psi->q,qrot2);
            psi->q = qrot2;
            rotmat = getrotmatrix(psi->q);
            matrix_x_vector(rotmat,sys.site[0],u);
            r2 = vector_inp(u,tetrasites[ipart][0]);
            if(r2>0.9999999) {
                foundquat=1;
                
            }
            count++;
            if(count==10000000) {
                foundquat=1;
            }
        } while(foundquat==0);
    }
    for(ipart=0; ipart<sys.npart; ipart++) {
        qprint(psi->q);
    }

    angle=PI*0./180.;
    vector_minus(psl->pts[0].r,psl->pts[1].r,rotax);
    r2=vector_inp(rotax,rotax);
    scalar_divide(rotax,sqrt(r2),rotax);
    qrot.q0 = cos(angle/2.0);
    qrot.q1 = rotax.x*sin(angle/2.0);
    qrot.q2 = rotax.y*sin(angle/2.0);
    qrot.q3 = rotax.z*sin(angle/2.0);
    quat_times(qrot,psl->pts[1].q,qrot2);
    psl->pts[1].q=qrot2;




    printf("Total energy system:            %lf\n",potential_energy(psl));
    printf("Total number of bonds:          %d\n",psl->nbonds);
    printf("\n");


    return;
}



void setup_octahedron(Slice *psl) {
    
    int ix,iy,iz,sign;
    int ipart,jpart,ipatch,isite;
    double phi,ru,scale;
    Pts *psi,*psj;
    double r2,l,sqrt2;
    vector dr;
    vector octasites[sys.npart][sys.nsites],N,M,u0,u1;


    printf("Setting all positions for the Octahedron and scaling bond lengths to sigma_LJ\n");

    ipart=0;
    for(ix=-1; ix<=1; ix+=2) {
        psl->pts[ipart].r.x = ix;
        psl->pts[ipart].r.z = 0;
        psl->pts[ipart].r.y = 0.;
        ipart++;
    }

    for(iy=-1; iy<=1; iy+=2) {
        psl->pts[ipart].r.y = iy;
        psl->pts[ipart].r.z = 0;
        psl->pts[ipart].r.x = 0.;
        ipart++;
    }

    for(iz=-1; iz<=1; iz+=2) {
        psl->pts[ipart].r.z = iz;
        psl->pts[ipart].r.y = 0;
        psl->pts[ipart].r.x = 0.;
        ipart++;
    }



    ru = 0.5*sqrt(2.);
    scale = sys.sigmaLJ*ru;
    for(ipart=0; ipart<sys.npart; ipart++) {
        scalar_times(psl->pts[ipart].r,scale,psl->pts[ipart].r);
    }


    printf("Positions of all particles\n");
    for(ipart=0; ipart<sys.npart; ipart++) {
        printf("%lf %lf %lf\n",psl->pts[ipart].r.x,psl->pts[ipart].r.y,psl->pts[ipart].r.z);
    }



    printf("Obtaining patches for all particles\n");
    for(ipart=0; ipart<sys.npart; ipart++) {
        psi = &(psl->pts[ipart]);
        ipatch =0;
        for(jpart=0; jpart<sys.npart; jpart++) {
            if(jpart!=ipart) {
                psj = &(psl->pts[jpart]);
                vector_minus(psj->r,psi->r,dr);
                r2 = vector_inp(dr,dr);
                l = sqrt(r2);
                if(l<(sys.sigmaLJ+0.001)) {
                    scalar_divide(dr,l,dr);
                    octasites[ipart][ipatch] = dr;
                    //r2=vector_inp(octasites[ipart][ipatch],octasites[ipart][ipatch]);
                    //printf("length patch %d particle %d: %lf\n",ipatch,ipart,sqrt(r2));
                    ipatch++;
                }
            }
        }
        //dprint(ipatch);
    }

    //find proper quaternion for every particle
    //  first rotate vector (0,0,1) in proper direction of particle to centre of polyhedron
    //  subsequently rotate one of the patches along rotate (0,0,1) axis, other should immediately follow

    vector rotax;
    tensor rotmat;
    int foundquat,count;
    double angle;
    quaternion qrot,qrot2;
    vector u,zvec,zvecipart;

    zvec.x=0; zvec.y=0; zvec.z=1.0;

    printf("Obtaining correct quaternion for rotation\n");
    for(ipart=0; ipart<sys.npart; ipart++) {
        //printf("Particle %d\n",ipart);
        psi = &(psl->pts[ipart]);
        scalar_divide(psi->r,-sqrt(vector_inp(psi->r,psi->r)),zvecipart);
        qrot = quatVecToVec(zvec,zvecipart);
        rotax = zvecipart;
        psi->q = qrot;
        angle = 0.001;
        qrot.q0 = cos(angle/2.0);
        qrot.q1 = rotax.x*sin(angle/2.0);
        qrot.q2 = rotax.y*sin(angle/2.0);
        qrot.q3 = rotax.z*sin(angle/2.0);
        count=0;
        do {
            foundquat=0;
            quat_times(qrot,psi->q,qrot2);
            psi->q = qrot2;
            rotmat = getrotmatrix(psi->q);
            matrix_x_vector(rotmat,sys.site[0],u);
            r2 = vector_inp(u,octasites[ipart][0]);
            if(r2>0.9999999) {
                foundquat=1;
                
            }
            count++;
            if(count==100000) {
                foundquat=1;
            }
        } while(foundquat==0);
    }

    angle=PI*0./180.;
    vector_minus(psl->pts[0].r,psl->pts[2].r,rotax);
    r2=vector_inp(rotax,rotax);
    scalar_divide(rotax,sqrt(r2),rotax);
    qrot.q0 = cos(angle/2.0);
    qrot.q1 = rotax.x*sin(angle/2.0);
    qrot.q2 = rotax.y*sin(angle/2.0);
    qrot.q3 = rotax.z*sin(angle/2.0);
    quat_times(qrot,psl->pts[2].q,qrot2);
    psl->pts[2].q=qrot2;


    printf("Total energy system:            %lf\n",potential_energy(psl));
    printf("Total number of bonds:          %d\n",psl->nbonds);
    printf("\n");


    return;
}


void setup_cube(Slice *psl) {
    
    int ix,iy,iz,sign;
    int ipart,jpart,ipatch,isite;
    double phi,ru,scale;
    Pts *psi,*psj;
    double r2,l,sqrt2;
    vector dr;
    vector cubesites[sys.npart][sys.nsites],N,M,u0,u1;


    printf("Setting all positions for the Octahedron and scaling bond lengths to sigma_LJ\n");

    ipart=0;
    for(ix=-1; ix<=1; ix+=2) {
        for(iy=-1; iy<=1; iy+=2) {
            for(iz=-1; iz<=1; iz+=2) {
                psl->pts[ipart].r.x = ix;
                psl->pts[ipart].r.z = iy;
                psl->pts[ipart].r.y = iz;
                ipart++;
            }
        }
    }


    ru = sqrt(0.25+0.25+0.25);
    scale = sys.sigmaLJ*ru/sqrt(1.+1.+1.);
    for(ipart=0; ipart<sys.npart; ipart++) {
        scalar_times(psl->pts[ipart].r,scale,psl->pts[ipart].r);
    }


    printf("Positions of all particles\n");
    for(ipart=0; ipart<sys.npart; ipart++) {
        printf("%lf %lf %lf\n",psl->pts[ipart].r.x,psl->pts[ipart].r.y,psl->pts[ipart].r.z);
    }



    printf("Obtaining patches for all particles\n");
    for(ipart=0; ipart<sys.npart; ipart++) {
        psi = &(psl->pts[ipart]);
        ipatch =0;
        for(jpart=0; jpart<sys.npart; jpart++) {
            if(jpart!=ipart) {
                psj = &(psl->pts[jpart]);
                vector_minus(psj->r,psi->r,dr);
                r2 = vector_inp(dr,dr);
                l = sqrt(r2);
                if(l<(sys.sigmaLJ+0.001)) {
                    scalar_divide(dr,l,dr);
                    cubesites[ipart][ipatch] = dr;
                    //r2=vector_inp(cubesites[ipart][ipatch],cubesites[ipart][ipatch]);
                    //printf("length patch %d particle %d: %lf\n",ipatch,ipart,sqrt(r2));
                    ipatch++;
                }
            }
        }
        //dprint(ipatch);
    }

    //find proper quaternion for every particle
    //  first rotate vector (0,0,1) in proper direction of particle to centre of polyhedron
    //  subsequently rotate one of the patches along rotate (0,0,1) axis, other should immediately follow

    vector rotax;
    tensor rotmat;
    int foundquat,count;
    double angle;
    quaternion qrot,qrot2;
    vector u,zvec,zvecipart;

    zvec.x=0; zvec.y=0; zvec.z=1.0;

    printf("Obtaining correct quaternion for rotation\n");
    for(ipart=0; ipart<sys.npart; ipart++) {
        //printf("Particle %d\n",ipart);
        psi = &(psl->pts[ipart]);
        scalar_divide(psi->r,-sqrt(vector_inp(psi->r,psi->r)),zvecipart);
        qrot = quatVecToVec(zvec,zvecipart);
        rotax = zvecipart;
        psi->q = qrot;
        angle = 0.001;
        qrot.q0 = cos(angle/2.0);
        qrot.q1 = rotax.x*sin(angle/2.0);
        qrot.q2 = rotax.y*sin(angle/2.0);
        qrot.q3 = rotax.z*sin(angle/2.0);
        count=0;
        do {
            foundquat=0;
            quat_times(qrot,psi->q,qrot2);
            psi->q = qrot2;
            rotmat = getrotmatrix(psi->q);
            matrix_x_vector(rotmat,sys.site[0],u);
            r2 = vector_inp(u,cubesites[ipart][0]);
            if(r2>0.9999999) {
                foundquat=1;
                
            }
            count++;
            if(count==100000) {
                foundquat=1;
            }
        } while(foundquat==0);
    }

    angle=PI*0./180.;
    vector_minus(psl->pts[0].r,psl->pts[2].r,rotax);
    r2=vector_inp(rotax,rotax);
    scalar_divide(rotax,sqrt(r2),rotax);
    qrot.q0 = cos(angle/2.0);
    qrot.q1 = rotax.x*sin(angle/2.0);
    qrot.q2 = rotax.y*sin(angle/2.0);
    qrot.q3 = rotax.z*sin(angle/2.0);
    quat_times(qrot,psl->pts[2].q,qrot2);
    psl->pts[2].q=qrot2;


    printf("Total energy system:            %lf\n",potential_energy(psl));
    printf("Total number of bonds:          %d\n",psl->nbonds);
    printf("\n");


    return;
}





void setup_icosahedron(Slice *psl) {
    
    int ix,iy,iz,sign;
    int ipart,jpart,ipatch,isite;
    double phi,ru,scale;
    Pts *psi,*psj;
    double r2,l;
    vector dr;
    vector icosites[sys.npart][sys.nsites],N,M,u0,u1;

    phi = (1.+sqrt(5.))/2.;

    printf("Setting all positions for the Icosahedron and scaling bond lengths to sigma_LJ\n");
    ipart=0;
    for(iy=-1; iy<=1; iy+=2) {
        for(iz=-1; iz<=1; iz+=2) {
            psl->pts[ipart].r.y = (1.)*(double)iy;
            psl->pts[ipart].r.z = (phi)*(double)iz;
            psl->pts[ipart].r.x = 0.;
            ipart++;
        }
    }

    for(ix=-1; ix<=1; ix+=2) {
        for(iy=-1; iy<=1; iy+=2) {
            psl->pts[ipart].r.x = (1.)*(double)ix;
            psl->pts[ipart].r.y = (phi)*(double)iy;
            psl->pts[ipart].r.z = 0.;
            ipart++;
        }
    }

    for(ix=-1; ix<=1; ix+=2) {
        for(iz=-1; iz<=1; iz+=2) {
            psl->pts[ipart].r.x = (phi)*(double)ix;
            psl->pts[ipart].r.z = (1.)*(double)iz;
            psl->pts[ipart].r.y = 0.;
            ipart++;
        }
    }


    ru = sqrt(sqrt(5.)*phi)*0.5;
    scale = sys.sigmaLJ*ru/sqrt(1.+(phi*phi));
    for(ipart=0; ipart<sys.npart; ipart++) {
        scalar_times(psl->pts[ipart].r,scale,psl->pts[ipart].r);
    }


    printf("Positions of all particles\n");
    for(ipart=0; ipart<sys.npart; ipart++) {
        printf("%lf %lf %lf\n",psl->pts[ipart].r.x,psl->pts[ipart].r.y,psl->pts[ipart].r.z);
    }



    printf("Obtaining correct quaternion for rotation\n");
    for(ipart=0; ipart<sys.npart; ipart++) {
        psi = &(psl->pts[ipart]);
        ipatch =0;
        for(jpart=0; jpart<sys.npart; jpart++) {
            if(jpart!=ipart) {
                psj = &(psl->pts[jpart]);
                vector_minus(psj->r,psi->r,dr);
                r2 = vector_inp(dr,dr);
                l = sqrt(r2);
                if(l<(sys.sigmaLJ+0.001)) {
                    scalar_divide(dr,l,dr);
                    icosites[ipart][ipatch] = dr;
                    //r2=vector_inp(icosites[ipart][ipatch],icosites[ipart][ipatch]);
                    //printf("length patch %d particle %d: %lf\n",ipatch,ipart,sqrt(r2));
                    ipatch++;
                }
            }
        }
        //dprint(ipatch);
    }

    //find proper quaternion for every particle
    //  first rotate vector (0,0,1) in proper direction of particle to centre of polyhedron
    //  subsequently rotate one of the patches along rotate (0,0,1) axis, other should immediately follow

    vector rotax;
    tensor rotmat;
    int foundquat,count;
    double angle;
    quaternion qrot,qrot2;
    vector u,zvec,zvecipart;

    zvec.x=0; zvec.y=0; zvec.z=1.0;

    for(ipart=0; ipart<sys.npart; ipart++) {
        //printf("Particle %d\n",ipart);
        psi = &(psl->pts[ipart]);
        scalar_divide(psi->r,-sqrt(vector_inp(psi->r,psi->r)),zvecipart);
        qrot = quatVecToVec(zvec,zvecipart);
        rotax = zvecipart;
        psi->q = qrot;
        angle = 0.001;
        qrot.q0 = cos(angle/2.0);
        qrot.q1 = rotax.x*sin(angle/2.0);
        qrot.q2 = rotax.y*sin(angle/2.0);
        qrot.q3 = rotax.z*sin(angle/2.0);
        count=0;
        do {
            foundquat=0;
            quat_times(qrot,psi->q,qrot2);
            psi->q = qrot2;
            rotmat = getrotmatrix(psi->q);
            matrix_x_vector(rotmat,sys.site[0],u);
            r2 = vector_inp(u,icosites[ipart][0]);
            if(r2>0.9999999) {
                foundquat=1;
                
            }
            count++;
            if(count==10000000) {
                foundquat=1;
            }
        } while(foundquat==0);
    }

    angle=PI*0.0/180.;
    vector_minus(psl->pts[0].r,psl->pts[8].r,rotax);
    r2=vector_inp(rotax,rotax);
    scalar_divide(rotax,sqrt(r2),rotax);
    qrot.q0 = cos(angle/2.0);
    qrot.q1 = rotax.x*sin(angle/2.0);
    qrot.q2 = rotax.y*sin(angle/2.0);
    qrot.q3 = rotax.z*sin(angle/2.0);
    quat_times(qrot,psl->pts[8].q,qrot2);
    psl->pts[8].q=qrot2;



    printf("Total energy system: %lf\n",potential_energy(psl));
    printf("Total number of bonds: %d\n",psl->nbonds);


    return;
}



void setup_dodecahedron(Slice *psl) {
    
    int ix,iy,iz,sign;
    int ipart,jpart,ipatch,isite;
    double phi,ru,scale;
    Pts *psi,*psj;
    double r2,l;
    vector dr;
    vector dodecasites[sys.npart][sys.nsites],N,M,u0,u1;

    phi = (sqrt(5.)-1.0)/2.;

    printf("Setting all positions for the Dodecahedron and scaling bond lengths to sigma_LJ\n");
    ipart=0;

    for(ix=-1; ix<=1; ix+=2) {
        for(iy=-1; iy<=1; iy+=2) {
            for(iz=-1; iz<=1; iz+=2) {
                psl->pts[ipart].r.x = (double)ix;
                psl->pts[ipart].r.y = (double)iy;
                psl->pts[ipart].r.z = (double)iz;
                ipart++;
            }
        }
    }


    for(iy=-1; iy<=1; iy+=2) {
        for(iz=-1; iz<=1; iz+=2) {
            psl->pts[ipart].r.x = 0.;
            psl->pts[ipart].r.y = (1.+phi)*(double)iy;
            psl->pts[ipart].r.z = (1.-(phi*phi))*(double)iz;
            ipart++;
        }
    }

    for(ix=-1; ix<=1; ix+=2) {
        for(iy=-1; iy<=1; iy+=2) {
            psl->pts[ipart].r.x = (1.+phi)*(double)ix;
            psl->pts[ipart].r.y = (1.-(phi*phi))*(double)iy;
            psl->pts[ipart].r.z = 0.;
            ipart++;
        }
    }

    for(ix=-1; ix<=1; ix+=2) {
        for(iz=-1; iz<=1; iz+=2) {
            psl->pts[ipart].r.x = (1.-(phi*phi))*(double)ix;
            psl->pts[ipart].r.y = 0.;
            psl->pts[ipart].r.z = (1.+phi)*(double)iz;
            ipart++;
        }
    }

    dprint(ipart);

    ru = sqrt(3.)*(sqrt(5.)+1)*0.25;
    scale = sys.sigmaLJ*ru/sqrt(3.);
    for(ipart=0; ipart<sys.npart; ipart++) {
        scalar_times(psl->pts[ipart].r,scale,psl->pts[ipart].r);
    }


    printf("Positions of all particles\n");
    for(ipart=0; ipart<sys.npart; ipart++) {
        printf("%lf %lf %lf\n",psl->pts[ipart].r.x,psl->pts[ipart].r.y,psl->pts[ipart].r.z);
    }



    printf("Obtaining correct quaternion for rotation\n");
    for(ipart=0; ipart<sys.npart; ipart++) {
        psi = &(psl->pts[ipart]);
        ipatch =0;
        for(jpart=0; jpart<sys.npart; jpart++) {
            if(jpart!=ipart) {
                psj = &(psl->pts[jpart]);
                vector_minus(psj->r,psi->r,dr);
                r2 = vector_inp(dr,dr);
                l = sqrt(r2);
                if(l<(sys.sigmaLJ+0.1)) {
                    scalar_divide(dr,l,dr);
                    dodecasites[ipart][ipatch] = dr;
                    //r2=vector_inp(dodecasites[ipart][ipatch],dodecasites[ipart][ipatch]);
                    //printf("length patch %d particle %d: %lf\n",ipatch,ipart,sqrt(r2));
                    ipatch++;
                }
            }
        }
        //dprint(ipatch);
    }

    //find proper quaternion for every particle
    //  first rotate vector (0,0,1) in proper direction of particle to centre of polyhedron
    //  subsequently rotate one of the patches along rotate (0,0,1) axis, other should immediately follow

    vector rotax;
    tensor rotmat;
    int foundquat,count;
    double angle;
    quaternion qrot,qrot2;
    vector u,zvec,zvecipart,utemp[3];

    u=dodecasites[0][0];
    vector_add(u,dodecasites[0][1],u);
    vector_add(u,dodecasites[0][2],u);
    l=sqrt(vector_inp(u,u));
    scalar_divide(u,l,u);

    zvec.x=0; zvec.y=0; zvec.z=1.0;
    qrot = quatVecToVec(u,zvec);
    rotmat=getrotmatrix(qrot);
    matrix_x_vector(rotmat,dodecasites[0][0],u);
    printf("ux %lf uy %lf uz %.16lf\n",u.x,u.y,u.z);

    for(ipart=0; ipart<sys.npart; ipart++) {
        //printf("Particle %d\n",ipart);
        psi = &(psl->pts[ipart]);
        l=-sqrt(vector_inp(psi->r,psi->r));
        scalar_divide(psi->r,l,zvecipart);
        qrot = quatVecToVec(zvec,zvecipart);
        rotax = zvecipart;
        psi->q = qrot;
        angle = 0.001;
        qrot.q0 = cos(angle/2.0);
        qrot.q1 = rotax.x*sin(angle/2.0);
        qrot.q2 = rotax.y*sin(angle/2.0);
        qrot.q3 = rotax.z*sin(angle/2.0);
        count=0;
        do {
            foundquat=0;
            quat_times(qrot,psi->q,qrot2);
            psi->q = qrot2;
            rotmat = getrotmatrix(psi->q);
            matrix_x_vector(rotmat,sys.site[0],u);
            r2 = vector_inp(u,dodecasites[ipart][0]);
            if(r2>0.999999) {
                foundquat=1;
                
            }
            count++;
            if(count==100000000) {
                foundquat=1;
            }
        } while(foundquat==0);
    }

    angle=PI*0.0/180.;
    vector_minus(psl->pts[0].r,psl->pts[8].r,rotax);
    r2=vector_inp(rotax,rotax);
    scalar_divide(rotax,sqrt(r2),rotax);
    qrot.q0 = cos(angle/2.0);
    qrot.q1 = rotax.x*sin(angle/2.0);
    qrot.q2 = rotax.y*sin(angle/2.0);
    qrot.q3 = rotax.z*sin(angle/2.0);
    quat_times(qrot,psl->pts[8].q,qrot2);
    psl->pts[8].q=qrot2;


    printf("Total energy system: %lf\n",potential_energy(psl));
    printf("Total number of bonds: %d\n",psl->nbonds);


    return;
}


void setup_truncatedicosahedron(Slice *psl) {
    
    int ix,iy,iz,sign;
    int ipart,jpart,ipatch,isite,jsite;
    double phi,ru,scale;
    Pts *psi,*psj;
    double r2,l;
    double C0,C1,C2,C3,C4;
    vector dr;
    vector truncicosasites[sys.npart][sys.nsites],N,M,u0,u1;
    vector rotax;
    tensor rotmat;
    int foundquat,count;
    double angle;
    quaternion qrot,qrot2;
    vector u,zvec,zvecipart;


    phi=(1.+sqrt(5.))/2.;

    printf("Setting all positions for the Truncated Icosahedron and scaling bond lengths to sigma_LJ\n");

    ipart=0;
    //x axis permutation
    for(ix=-1; ix<=1; ix+=2) {
        for(iz=-1; iz<=1; iz+=2) {
            psl->pts[ipart].r.x = ix*3.*phi;
            psl->pts[ipart].r.y = 0;
            psl->pts[ipart].r.z = iz;
            ipart++;
        }
    }
    for(ix=-1; ix<=1; ix+=2) {
        for(iy=-1; iy<=1; iy+=2) {
            for(iz=-1; iz<=1; iz+=2) {
                psl->pts[ipart].r.x = ix*(1.+2.*phi);
                psl->pts[ipart].r.y = iy*phi;
                psl->pts[ipart].r.z = iz*2.;
                ipart++;
            }
        }
    }
    for(ix=-1; ix<=1; ix+=2) {
        for(iy=-1; iy<=1; iy+=2) {
            for(iz=-1; iz<=1; iz+=2) {
                psl->pts[ipart].r.x = ix*(2.+phi);
                psl->pts[ipart].r.y = iy*2.*phi;
                psl->pts[ipart].r.z = iz;
                ipart++;
            }
        }
    }


    //y axis permutation
    for(ix=-1; ix<=1; ix+=2) {
        for(iy=-1; iy<=1; iy+=2) {
            psl->pts[ipart].r.x = ix;
            psl->pts[ipart].r.y = iy*3.*phi;
            psl->pts[ipart].r.z = 0;
            ipart++;
        }
    }
    for(ix=-1; ix<=1; ix+=2) {
        for(iy=-1; iy<=1; iy+=2) {
            for(iz=-1; iz<=1; iz+=2) {
                psl->pts[ipart].r.x = ix*2.;
                psl->pts[ipart].r.y = iy*(1.+2.*phi);
                psl->pts[ipart].r.z = iz*phi;
                ipart++;
            }
        }
    }
    for(ix=-1; ix<=1; ix+=2) {
        for(iy=-1; iy<=1; iy+=2) {
            for(iz=-1; iz<=1; iz+=2) {
                psl->pts[ipart].r.x = ix;
                psl->pts[ipart].r.y = iy*(2.+phi);
                psl->pts[ipart].r.z = iz*2.*phi;
                ipart++;
            }
        }
    }

    //z axis permutation
    for(iy=-1; iy<=1; iy+=2) {
        for(iz=-1; iz<=1; iz+=2) {
            psl->pts[ipart].r.x = 0;
            psl->pts[ipart].r.y = iy;
            psl->pts[ipart].r.z = iz*3.*phi;
            ipart++;
        }
    }
    for(ix=-1; ix<=1; ix+=2) {
        for(iy=-1; iy<=1; iy+=2) {
            for(iz=-1; iz<=1; iz+=2) {
                psl->pts[ipart].r.x = ix*phi;
                psl->pts[ipart].r.y = iy*2.;
                psl->pts[ipart].r.z = iz*(1.+2.*phi);
                ipart++;
            }
        }
    }
    for(ix=-1; ix<=1; ix+=2) {
        for(iy=-1; iy<=1; iy+=2) {
            for(iz=-1; iz<=1; iz+=2) {
                psl->pts[ipart].r.x = ix*2.*phi;
                psl->pts[ipart].r.y = iy;
                psl->pts[ipart].r.z = iz*(2.+phi);
                ipart++;
            }
        }
    }


    dprint(ipart);

    ru = sqrt(0.5*0.5 + (3.*(1.+sqrt(5.))/4.)*(3.*(1.+sqrt(5.))/4.));
    scale = sys.sigmaLJ*ru/sqrt((1.+(3.*phi)*(3.*phi)));
    for(ipart=0; ipart<sys.npart; ipart++) {
        scalar_times(psl->pts[ipart].r,scale,psl->pts[ipart].r);
    }


    //printf("Positions of all particles\n");
    //for(ipart=0; ipart<sys.npart; ipart++) {
    //    printf("%d %lf %lf %lf\n",ipart,psl->pts[ipart].r.x,psl->pts[ipart].r.y,psl->pts[ipart].r.z);
    //}


    printf("Obtaining correct quaternion for rotation\n");
    for(ipart=0; ipart<sys.npart; ipart++) {
        psi = &(psl->pts[ipart]);
        ipatch =0;
        //dprint(ipart);
        for(jpart=0; jpart<sys.npart; jpart++) {
            if(jpart!=ipart) {
                psj = &(psl->pts[jpart]);
                vector_minus(psj->r,psi->r,dr);
                r2 = vector_inp(dr,dr);
                l = sqrt(r2);
                if(l<(sys.sigmaLJ+0.1)) {
                    //if(ipart==0) {
                    //    printf("ipart %d jpart %d are connected\n",ipart,jpart);
                    //}
                    scalar_divide(dr,l,dr);
                    truncicosasites[ipart][ipatch] = dr;
                    ipatch++;
                }
            }
        }
        //printf("ipart %d has %d patches\n",ipart,ipatch);
    }



    zvec.x=0.;
    zvec.y=0.;
    zvec.z=1.;

    int correctsite,countcorrectsite;
    for(ipart=0; ipart<sys.npart; ipart++) {
        psi = &(psl->pts[ipart]);
        l=-sqrt(vector_inp(psi->r,psi->r));
        scalar_divide(psi->r,l,zvecipart);
        qrot = quatVecToVec(zvec,zvecipart);
        rotax = zvecipart;
        psi->q = qrot;
        angle = 2.*PI/10000000.;
        qrot.q0 = cos(angle/2.0);
        qrot.q1 = rotax.x*sin(angle/2.0);
        qrot.q2 = rotax.y*sin(angle/2.0);
        qrot.q3 = rotax.z*sin(angle/2.0);
        count=0;
        for(isite=0; isite<sys.nsites; isite++) {
            countcorrectsite=0;
            for(jsite=0; jsite<sys.nsites; jsite++) {
                if(isite!=jsite) {
                    r2=vector_inp(truncicosasites[ipart][isite],truncicosasites[ipart][jsite]);
                    if(((180./PI)*acos(r2))>110.) {
                        countcorrectsite++;
                    }
                }
            }
            if(countcorrectsite==2) {
                correctsite=isite;
                //printf("ipart %d correctsite %d\n",ipart,isite);
                break;
            }
        }

        do {
            foundquat=0;
            quat_times(qrot,psi->q,qrot2);
            psi->q = qrot2;
            rotmat = getrotmatrix(psi->q);
            matrix_x_vector(rotmat,sys.site[0],u);
            r2 = vector_inp(u,truncicosasites[ipart][correctsite]);
            if(r2>0.9999999) {
                foundquat=1;
                //printf("found by reaching small enough r2\n");
            }
            count++;
            if(count==100000000) {
                //printf("not found by reaching small enough r2\n");
                foundquat=1;
            }
        } while(foundquat==0);
    }

    

    angle=PI*0.0/180.;
    vector_minus(psl->pts[0].r,psl->pts[1].r,rotax);
    r2=vector_inp(rotax,rotax);
    scalar_divide(rotax,sqrt(r2),rotax);
    qrot.q0 = cos(angle/2.0);
    qrot.q1 = rotax.x*sin(angle/2.0);
    qrot.q2 = rotax.y*sin(angle/2.0);
    qrot.q3 = rotax.z*sin(angle/2.0);
    quat_times(qrot,psl->pts[1].q,qrot2);
    psl->pts[1].q=qrot2;




    return;
}



void print_tis_setup() {

    int istate, totalrep,type;

    printf("TIS setup:\n");

    totalrep = 0;
    printf("Total number of states defined: %d\n", path.nstates);
    printf("Total number of replicas defined: ");
    for(istate=0; istate<path.nstates; istate++) {
        printf("%d + ",state[istate].nrep);
        totalrep+=state[istate].nrep;
    }
    printf("= %d\n",totalrep);

    printf("Current initial state:                         %d\n",path.initial_state);
    printf("Current initial state according to in_state(): %d\n", in_state(&slice[0]));
    printf("Current middle state according to in_state():  %d\n", in_state(&slice[(int)((path.nslices-1)/2.0)]));
    printf("Current final state according to in_state():   %d\n", in_state(&slice[path.nslices-1]));
    printf("Current replica: %d\n", path.current_replica);
    if(path.current_replica==0) {
        type =analyse_state_i(slice,replica[path.current_replica],path.nslices, path.initial_state);
    }
    else {
        type =analyse(slice,replica[path.current_replica],path.nslices, path.initial_state);
    }
    printf("Type of current path: %d\n", type);
    printf("Length current path: %d\n", path.nslices);


    printf("\n");

    return;
}





void stats_init() {

    int istate,irep,iaver,iacc;

    printf("Initializing stats\n");
    for(istate=0; istate<MAXSTATES; istate++) {
        for(irep=0; irep<MAXREPLICA; irep++) {
            for(iaver=0; iaver<NSTAT; iaver++) {
                path.block_stats[istate][irep].aver[iaver].now=0;
                path.block_stats[istate][irep].aver[iaver].sum=0;
                path.block_stats[istate][irep].aver[iaver].sumsq=0;
                path.block_stats[istate][irep].aver[iaver].n=0;
                path.final_stats[istate][irep].aver[iaver].now=0;
                path.final_stats[istate][irep].aver[iaver].sum=0;
                path.final_stats[istate][irep].aver[iaver].sumsq=0;
                path.final_stats[istate][irep].aver[iaver].n=0;
            }
        }
    }

    for(istate=0; istate<MAXSTATES; istate++) {
        for(irep=0; irep<MAXREPLICA; irep++) {
            for(iacc=0; iacc<NACC; iacc++) {
                path.block_stats[istate][irep].mcacc[iacc].acc=0;
                path.block_stats[istate][irep].mcacc[iacc].tries=0;
                path.block_stats[istate][irep].mcacc[iacc].ratio=0;
                path.final_stats[istate][irep].mcacc[iacc].acc=0;
                path.final_stats[istate][irep].mcacc[iacc].tries=0;
                path.final_stats[istate][irep].mcacc[iacc].ratio=0;
            }
        }
    }

    for(istate=0; istate<MAXSTATES; istate++) {
        for(irep=0; irep<state[istate].nrep; irep++) {
            sprintf(path.block_stats[istate][irep].mcacc[0].name,"shots replica %3d ",irep);
            sprintf(path.block_stats[istate][irep].mcacc[1].name,"swaps replica %3d ",irep);
            sprintf(path.block_stats[istate][irep].mcacc[2].name,"swap  states  %3d ",irep);
            sprintf(path.block_stats[istate][irep].mcacc[3].name,"revs  replica %3d ",irep);
            sprintf(path.final_stats[istate][irep].mcacc[0].name,"shots replica %3d ",irep);
            sprintf(path.final_stats[istate][irep].mcacc[1].name,"swaps replica %3d ",irep);
            sprintf(path.final_stats[istate][irep].mcacc[2].name,"swap  states  %3d ",irep);
            sprintf(path.final_stats[istate][irep].mcacc[3].name,"revs  replica %3d ",irep);
        }
    }


    printf("Done initializing stats\n");

    return;
}







void traj_input() {

    char dum[40];
    char filename[100];
    FILE *fp;
    int i,j,type;
	Slice *psl;

	printf("\nIn traj_input()\n");
	printf("Reading in previous trajectory\n");
 
    sprintf(filename,"trajectory.inp");
    if ((fp = fopen(filename,"r"))==NULL){
        printf("input: can't open %s\n",filename);
        return;
    }
    else {
        fscanf(fp,"%d %s",&path.nslices,dum); 
        fscanf(fp,"%d %s",&path.current_replica,dum); 
        fscanf(fp,"%d %s",&path.initial_state,dum); 
        fscanf(fp,"%lf %lf %lf %s",&sys.boxl.x,&sys.boxl.y,&sys.boxl.z,dum);
        for(j=0;j<path.nslices;j++) {
            for(i=0;i<sys.npart;i++){
                fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",
                  &(slice[j].pts[i].r.x),&(slice[j].pts[i].r.y),&(slice[j].pts[i].r.z),
                  &(slice[j].pts[i].q.q0),&(slice[j].pts[i].q.q1),&(slice[j].pts[i].q.q2),&slice[j].pts[i].q.q3);
            }
        }
        fclose(fp);
    }  


    //initializing numbers not stored, but easily retraced
    path.current_gsen = state[path.initial_state-1].target.energy;
    path.nreplica = state[path.initial_state-1].nrep;
    path.final_state = in_state(&slice[path.nslices-1]);

    //if(path.initial_state==1) {
    //    state[0].target.energy = -slice[0].nbonds*sys.epsilonP;
    //}

    for(i=0;i<path.nslices;i++) {
        psl = &(slice[i]);
        create_all_rc(&slice[i]);
        //print_rc(&slice[i],path.initial_state);
    }


	path.current_gsen = state[path.initial_state-1].target.energy;
    printf("Warm start information:\n");
    printf("       number of slices:       %d\n",path.nslices);
    printf("       initial state:          %d\n",in_state(&slice[0]));
    printf("       middle state:           %d\n",in_state(&slice[(int)(path.nslices/2.0)]));
    printf("       final state:            %d\n",in_state(&slice[path.nslices-1]));
    printf("       current initial state:  %d\n",path.initial_state);
    printf("       current nreplica:       %d\n",path.nreplica);
    printf("       current replica:        %d\n",path.current_replica);
    printf("       current gsen:           %lf\n",path.current_gsen);
        

	for(i=0; i<path.nreplica; i++) {
		replica[i] = &state[path.initial_state-1].srep[i];
	}

    if(path.current_replica==0) {
        type =analyse_state_i(slice,replica[path.current_replica],path.nslices, path.initial_state);
        replica[path.current_replica]->type = type;
    }
    else {
        type =analyse(slice,replica[path.current_replica],path.nslices, path.initial_state);
        replica[path.current_replica]->type = type;
    }

    replica[path.current_replica]->pathlen=path.nslices;
  
	printf("Done with traj_input\n\n");

    return;
}



void traj_output() {

    FILE *fp;
    char filename[240];
    int i,j; 
    Replica *prep;

    prep = replica[path.current_replica];
  
    sprintf(filename,"trajectory.out");
    if ((fp = fopen(filename,"w"))==NULL) {
        printf("output: can't open %s\n",filename);
        return;
    }
    else {
        fprintf(fp,"%d slices\n",prep->pathlen); 
        fprintf(fp,"%d current_replica\n",path.current_replica); 
        fprintf(fp,"%d initial_state\n",path.initial_state); 
        fprintf(fp,"%.16lf %.16lf %.16lf boxl,dr\n",sys.boxl.x,sys.boxl.y,sys.boxl.z);
        for(j=0;j<path.nslices;j++) {
            for(i=0;i<sys.npart;i++){
                fprintf(fp,"%12.16lf %12.16lf %12.16lf %12.16lf %12.16lf %12.16lf %12.16lf\n",
                  slice[j].pts[i].r.x,slice[j].pts[i].r.y,slice[j].pts[i].r.z,
                  slice[j].pts[i].q.q0,slice[j].pts[i].q.q1,slice[j].pts[i].q.q2,slice[j].pts[i].q.q3);
            }
        }
        fclose(fp);
    }
	
    return;
}






void dos_input() {

    FILE *fp;
    char filename[30];
    int i,j,idum; 

    //Wang-Landau Dos
    sprintf(filename,"dos_all.inp");
    if ((fp = fopen(filename,"r"))==NULL){
        printf("output: can't open %s\n",filename);
        return;
    } 
    else {
        printf("reading dos data\n");
        for(i=0;i<path.nstates;i++) {
            fscanf(fp,"%d",&state[i].nrep);
            fscanf(fp,"%lf",&state[i].scalefactor);
            for(j=0;j<state[i].nrep;j++) {
                fscanf(fp,"%d %lf %lf\n",&idum,&state[i].srep[j].lambda,&state[i].srep[j].dos);
                printf("%d %g %g\n",j,state[i].srep[j].lambda,state[i].srep[j].dos);
            }
            fscanf(fp,"");
        }
        fclose(fp);
    }



	return;
}





void dos_output() {

	FILE *fp;
    int i,j;
	
    if ((fp = fopen("dos_all.dat","w"))==NULL){
        printf("output: can't open file.dat\n");
        return;
    }
    else {
        for(i=0;i<path.nstates;i++) {
            fprintf(fp,"%d\n",state[i].nrep);
            fprintf(fp,"%lf\n",state[i].scalefactor);
            for(j=0;j<state[i].nrep;j++) {
                fprintf(fp,"%d %.16lf %.16lf\n",j,state[i].srep[j].lambda,state[i].srep[j].dos);
            }
            fprintf(fp,"\n");
        }
        fclose(fp);
    }
	
    return;
}


void target_input() {

	FILE *fp;
    int istate,ipart;
	
    if ((fp = fopen("target.inp","r"))==NULL){
        printf("output: can't open target.inp\n");
        return;
    }
    else {
        fscanf(fp,"%d\n",&path.nstates);
        printf("%d number of states from target.out\n",path.nstates);
        for(istate=0;istate<path.nstates;istate++) {
            fscanf(fp,"%d %d %d %lf %lf\n",
                    &state[istate].target.nbonds,
                    &state[istate].target.maxclustersize,
                    &state[istate].target.misalignments,
                    &state[istate].target.energy,
                    &state[istate].volume_op);
            print_rc(&state[istate].target,istate+1);
            for(ipart=0;ipart<sys.npart;ipart++) {
                fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf\n",
                        &state[istate].target.pts[ipart].r.x,
                        &state[istate].target.pts[ipart].r.y,
                        &state[istate].target.pts[ipart].r.z,
                        &state[istate].target.pts[ipart].q.q0,
                        &state[istate].target.pts[ipart].q.q1,
                        &state[istate].target.pts[ipart].q.q2,
                        &state[istate].target.pts[ipart].q.q3);
                //vprint(state[istate].target.pts[ipart].r);
                //qprint(state[istate].target.pts[ipart].q);
            }
            fscanf(fp,"\n");
        }
        fclose(fp);
    }


    return;
}



void target_output() {

	FILE *fp;
    int istate,ipart;
	
    if ((fp = fopen("target.out","w"))==NULL){
        printf("output: can't open target.out\n");
        return;
    }
    else {
        fprintf(fp,"%d\n",path.nstates);
        for(istate=0;istate<path.nstates;istate++) {
            fprintf(fp,"%d %d %d %lf %lf\n",
                    state[istate].target.nbonds,
                    state[istate].target.maxclustersize,
                    state[istate].target.misalignments,
                    state[istate].target.energy,
                    state[istate].volume_op);
            for(ipart=0;ipart<sys.npart;ipart++) {
                fprintf(fp,"%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf\n",
                        state[istate].target.pts[ipart].r.x,
                        state[istate].target.pts[ipart].r.y,
                        state[istate].target.pts[ipart].r.z,
                        state[istate].target.pts[ipart].q.q0,
                        state[istate].target.pts[ipart].q.q1,
                        state[istate].target.pts[ipart].q.q2,
                        state[istate].target.pts[ipart].q.q3);
            }
            fprintf(fp,"\n");
        }
        fclose(fp);
    }

    return;
}



void warm_start_gen() {

    target_output();
	traj_output();
	dos_output();

	return;
}





void conf_output(Slice *psl) {

    int ipart, isite;
    FILE *fp;


    if((fp=fopen("conf.out","w"))==NULL) {;
        printf("Warning: can not open conf.out\n");
    }
    else {
        fprintf(fp,"%d %d %lf\n", sys.npart, sys.nsites, sys.boxl.x);
        for(ipart=0; ipart<sys.npart; ipart++) {
            fprintf(fp,"%lf %lf %lf ", psl->pts[ipart].r.x, psl->pts[ipart].r.y, psl->pts[ipart].r.z);
            fprintf(fp,"%lf %lf %lf %lf\n", psl->pts[ipart].q.q0, psl->pts[ipart].q.q1, psl->pts[ipart].q.q2, psl->pts[ipart].q.q3);
        }
    }
    fclose(fp);

    return;
}




void conf_input(Slice *psl) {

    int ipart, isite,npart,nsites;
    double boxl;
    FILE *fp;


    if((fp=fopen("conf.inp","r"))==NULL) {;
        printf("Warning: can not open conf.out\n");
    }
    else {
        fscanf(fp,"%d %d %lf\n", &npart, &nsites, &boxl);
        for(ipart=0; ipart<npart; ipart++) {
            fscanf(fp,"%lf %lf %lf ", &psl->pts[ipart].r.x, &psl->pts[ipart].r.y, &psl->pts[ipart].r.z);
            fscanf(fp,"%lf %lf %lf %lf\n", &psl->pts[ipart].q.q0, &psl->pts[ipart].q.q1, &psl->pts[ipart].q.q2, &psl->pts[ipart].q.q3);
        }
    }
    fclose(fp);

    if(npart!=sys.npart) {
        printf("Warning: number of particles in system not same as in conf.inp\n");
    }
    if(nsites!=sys.nsites) {
        printf("Warning: number of sites in system not same as in conf.inp\n");
    }
    if(boxl!=sys.boxl.x) {
        printf("Warning: boxlength in system not same as in conf.inp\n");
    }


    return;
}


quaternion quatVecToVec(vector vec1, vector vec2) {

    quaternion qrot;
    vector rotvec,xvec,yvec;
    double angle,s,invs,l;

    angle=vector_inp(vec1,vec2);
    if(angle>=1) {
        qrot.q0=1;
        qrot.q1=0;
        qrot.q2=0;
        qrot.q3=0;
    }

    else if(angle < -0.999999) {
        //xvec.x=1;
        //xvec.y=0;
        //xvec.z=0;
        //vector_cross(xvec,vec2,rotvec);   
        //l=sqrt(vector_inp(rotvec,rotvec));
        //if(l==0) {
        //  yvec.x=0;
        //  yvec.y=1;
        //  yvec.z=0;
        //  vector_cross(yvec,vec2,rotvec);
        //  l=sqrt(vector_inp(rotvec,rotvec));
        //}
        //scalar_divide(rotvec,l,rotvec);
        //qrot.r=PI;
        //qrot.u=rotvec;
        //scdivide_quat(qrot,sqrt(quat_inp(qrot,qrot)),qrot);
        qrot.q0=0;
        qrot.q1=0;
        qrot.q2=1;
        qrot.q3=0;
    }

    else {
        s=sqrt(2.*(1.+angle));  
        invs=1./s;
        vector_cross(vec1,vec2,rotvec);
        qrot.q0=0.5*s;
        scalar_times(rotvec,invs,rotvec);
        qrot.q1=rotvec.x;
        qrot.q2=rotvec.y;
        qrot.q3=rotvec.z;
        scdivide_quat(qrot,sqrt(quat_inp(qrot,qrot)),qrot);
    }
                

    return qrot;
}

