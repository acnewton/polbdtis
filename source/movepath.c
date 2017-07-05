#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include "path.h"


void create_all_rc(Slice *psl) {

    int istate,ipart;

    psl->energy = potential_energy(psl);

    for(istate=0; istate<path.nstates; istate++) {
        if(istate==0) {
            psl->order_parameter[istate]= -psl->energy;
        }
        else if((istate+1)==path.initial_state) {
            psl->order_parameter[istate]=psl->enmaxclustop - state[istate].target.enmaxclust;
        }
        else {
            //psl->order_parameter[istate]=psl->enmaxclust - state[istate].target.enmaxclust;
            //psl->order_parameter[istate]=psl->energy - state[istate].target.energy;
            psl->order_parameter[istate]=psl->enmaxclust - state[istate].target.enmaxclust;
        }
    }

    return;
}
    


double get_rc(Slice *psl, int state_index) {

    return psl->order_parameter[state_index-1];

}



double print_rc(Slice *psl,int state_index) {

    double op;
    int ipart;
 
    printf("state %d instate %d en %lf enmaxclust %lf enmaxop %lf mindist %lf nbonds %d misalgns %d maxclust %d bondnotbreak %d op %lf\n",
            state_index,
            in_state(psl),
            psl->energy,
            psl->enmaxclust,
            psl->enmaxclustop,
            psl->mindist,
            psl->nbonds,
            psl->misalignments,
            psl->maxclustersize,
            psl->bonddidnotbreak,
            psl->order_parameter[state_index-1]);
    //for(ipart=0; ipart<sys.npart; ipart++) {
    //    printf("d%d %lf ",ipart,psl->opdist[ipart]);
    //}
    //printf("\n");


    return op;
}




int in_state(Slice *psl) {

    int istate;

    if(psl->nbonds==0) {
        if(psl->mindist>sys.rcutoffsq) {
            return 1;
        }
    }

    for(istate=1; istate<path.nstates; istate++) {
        //printf("State %d\n",istate);
        if(psl->mindist>sys.rcutoffsq) {
            //printf("mindist %lf\n",psl->mindist);
            if(psl->maxclustersize==state[istate].target.maxclustersize) {
                //printf("maxclustersize %d\n",psl->maxclustersize);
                if(psl->bonddidnotbreak==0) {
                    //printf("bondbreak\n",psl->bonddidnotbreak);
                    if(psl->nbonds==state[istate].target.nbonds) {
                        //printf("nbonds %d\n",psl->nbonds);
                        if(psl->misalignments==state[istate].target.misalignments) {
                            //printf("misalignments %d\n",psl->misalignments);
                            //if(psl->order_parameter[istate] < state[istate].volume_op) {
                            if((psl->enmaxclust-state[istate].target.enmaxclust) < state[istate].volume_op) {
                                //printf("op %lf\n",psl->order_parameter[istate]);
                                return istate+1;
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}




int trajectory(Slice *psl, int maxlength) {

    int i,istate;

    for (i=1;i<maxlength;i++) {
        psl[i] = psl[i-1]; 
        propagate_bd(&psl[i]);
        create_all_rc(&(psl[i]));
        if (in_state(&psl[i])>0) {
            return i;
        }
    }

    printf("maxlength %d reached ",maxlength);
    print_rc(&psl[maxlength-1],path.initial_state);

    return 0;
}



int in_upper_window(Slice *psl, Replica *prep, int state_index) {

    int in;
    double energy;

    if(prep->index>1) {
        energy =get_rc(psl,state_index);
        in = (energy>prep->lambda);
    }

    else {
        in = (in_state(psl)!=state_index);
    }
        
    return in; 

    //return (psl->order_parameter[state_index-1]>prep->lambda);
}


int trajectory_state_i(Slice *psl, Replica *prep, int maxlength, int state_index) {

    int i,istate;
    static int count,tries;

    for (i=1;i<maxlength;i++) {
        psl[i] = psl[i-1]; 
        propagate_bd(&psl[i]);
        create_all_rc(&(psl[i]));
        if (in_upper_window(&psl[i], prep,state_index)) {
            return i;
        }
    }

    return 0;
}



int analyse ( Slice *psl, Replica *prep, int len, int state_index) {

    int i,initial_state,final_state;
 
    //printf("In analyse()\n");
    //printf("Pathlength: %d\n",path.nslices);
    initial_state=in_state(&psl[0]);

    final_state=in_state(&psl[len-1]);

    if (initial_state==state_index) {
        if (final_state==state_index) {
            for (i=0;i<len;i++) {
                if (in_upper_window(&psl[i], prep,state_index)){
                    //printf("Out analyse()\n");
                    return 1;
                }
            }
        //printf("path did not reach lambda %g\n",prep->lambda);
        //print_rc(&psl[0],state_index);
        //print_rc(&psl[len-1],state_index);
            return 0;
        } 
        else {
            if (final_state==0) {
                printf("final state unknown\n");
                print_rc(&psl[0],state_index);
                print_rc(&psl[len-1],state_index);
                //printf("Out analyse()\n");
                return 0;
            }
            for (i=0;i<len;i++) {
                if (in_upper_window(&psl[i], prep,state_index)) {
                    //printf("intermediate slice %d is larger than lambda %g \n",i,prep->lambda);
                    //printf("Out analyse()\n");
                    return 2;
                }
            }
            //printf("warning: type 2 but path did not reach lambda\n");
            //dprint(initial_state);
            //dprint(final_state);
            //dprint(path.current_replica);
            //for(i=0; i<prep->pathlen; i++) {
            //  print_rc(&(trial[i]),state_index);  
            //  }
            //printf("Out analyse()\n");
            return 0;
        }
    }

    printf("analyse: path of length %d corrupted in replica %d, initial %d  final %d, state index %d \n",len, prep->index,initial_state,final_state,state_index);

    print_rc(&psl[0],state_index);
    print_rc(&psl[len-1],state_index);

    //printf("Out analyse()\n");
    return 0; 
}
  



int analyse_state_i(Slice *psl, Replica *prep, int len, int state_index) {

    int i,visitedA,visited2;
 
    visitedA=visited2 =0;
    if (in_upper_window(&psl[0], prep,state_index)){
        for (i=1;i<len;i++) {
            if (in_state(&psl[i])==state_index){
                visitedA =1;
            }
            if ((visitedA) && (in_upper_window(&psl[i], prep,state_index))) {
                visited2 =1;
            }
        }
        if ((visitedA) && (in_upper_window(&psl[len-1], prep,state_index))){
            //printf("path OK\n");
            return 1;
        }
        if ((visitedA) && (visited2) ) {
            printf ("path starts in 1, visits A and then 1 but doesn't end in 1\n");
            return 1;
        }
    }

    for(i=0;i<prep->pathlen;i++) {
        print_rc(&slice[i],state_index);
        printf("irep %d lambda %g energy %lf nbonds %d current_gsen %lf targetbonds %d",prep->index, prep->lambda, slice[i].energy, slice[i].nbonds, state[path.initial_state-1].target.energy, state[path.initial_state-1].target.nbonds );
        dprint(in_state(&slice[i]));
    }     
    printf("analyse i: path corrupted in replica %d   \n",prep->index);
    dprint(path.initial_state);

    return 0; 
}




int shoot_oneway(Replica *prep) {

    int i,j,islice,index,len,pathlen,maxlength,trial_pathlen,type;
    int initial_state,final_state;

    //printf("In shoot\n");

    if (prep->index==0) {
        //printf("Failed shoot because of prep index 0\n");
        return 0;
    }
  
    //do not want to keep minus interface statistics for shooting
    path.block_stats[path.initial_state-1][prep->index].mcacc[0].tries++;

    index =0;
    for (i=0;i<path.nslices;i++) {
        if (in_upper_window(&slice[i], prep, path.initial_state)) {
            index = i;
            i = prep->pathlen;
        }
    }
  
    //sanity check
    if (in_upper_window(&slice[index], prep, path.initial_state)==0) {
        printf("Failed shoot at sanity check\n");
        return 0;
    }
  
    for(i=0;i<=index;i++) {
        trial[i]=slice[i];
    }

    maxlength= MAXSLICES-index;
    len = trajectory(&trial[index],maxlength);
    if (len==0) {
        printf("Shoot oneway: trajectory too long. State %d\n",path.initial_state);
        if(sys.findnewstates==1) {
            printf("This trial trajectory might be ending up in a candidate state\n");
            find_minimum(&trial[maxlength-1]);
            j=in_state(&trial[maxlength-1]);
            if(j!=0) {
                printf("Via find_minimum() state %d is found\n",j);
                return 0;
            }
            if(add_newstate(&trial[maxlength-1])) {
                printf("New state found\n");
                //printf("Now checking if old path already visited new state\n");
                for(islice=0; islice<path.nslices; islice++) {
                    create_all_rc(&slice[islice]);
                }
                j=in_state(&slice[islice]);
                if( j!=path.initial_state && j!=0) {
                    printf("Found newly added state in old trajectory at slice %d. Breaking trajectory\n",
                            islice);
                    type = analyse(slice,prep,islice, path.initial_state);
                    if( type==0 ) {
                        printf("Warning: something went wrong with adding new state.\n");
                        dprint(path.initial_state);
                        dprint(initial_state);
                        dprint(final_state);
                        dprint(prep->index);
                        return 0;
                    }
                    prep->type =type;
                    path.nslices=islice+1;
                    prep->pathlen= islice+1;
                    final_state=in_state(&slice[path.nslices-1]);
                    path.final_state = final_state;
                    return 0;
                }
                return 0;
            }
            printf("No new state found\n");
            return 0;
        }
        return 0;
    }
    
    trial_pathlen = index+len+1;
    if(trial_pathlen<3) {
        //printf("path is too short: %d\n",trial_pathlen);
        printf("Failed shoot because path too short: %d\n",trial_pathlen);
        return 0;
    }

    type = analyse(trial,prep,trial_pathlen, path.initial_state);
    if( type==0 ) {
        printf("Shoot Oneway: wrong type after shooting.\n");
        dprint(path.initial_state);
        initial_state = in_state(&(trial[0]));
        dprint(trial[0].nbonds);
        gprint(trial[0].energy);
        dprint(initial_state);
        dprint(final_state);
        dprint(prep->index);

        //printf("Out shoot\n");
        return 0;
    }


    for(i=index+1;i<trial_pathlen;i++) {
        slice[i] = trial[i];
    }
  
    final_state=in_state(&slice[trial_pathlen-1]);
    path.final_state = final_state;

    prep->pathlen= trial_pathlen;
    path.nslices = trial_pathlen;
    prep->type =type;

    path.block_stats[path.initial_state-1][prep->index].mcacc[0].acc++;
    //printf("Out shoot\n");

    return 1;
}




int swap_replica(int irep, int jrep) {

    Replica *prepi,*prepj;
    Slice *swaptrial, h;
    int i,pathlen,ABexchOK,ireptype,jreptype,type;
    double  aux;
    //printf("In swaprep\n");

    prepi = replica[irep];
    path.block_stats[path.initial_state-1][prepi->index].mcacc[1].tries++;
    prepj = replica[jrep];

    //if(path.fixedbias==0) {
        prepi->dos += state[path.initial_state-1].scalefactor;
    //}
  
    if (jrep<0) {
        //printf("Out swaprep\n");
        return 0;
    }

    if ( (irep==0) || (jrep==0) ) {
        //printf("Out swaprep\n");
        return swap_replica_0(irep, jrep);
    }

    if (( irep==path.nreplica-1 ) && ( jrep==path.nreplica )) {
        //printf("Out swaprep\n");
        return 0;
    }
    
    if ( jrep>path.nreplica-1 ) {
        //printf("Out swaprep\n");
        return 0;
    }

    aux = prepi->dos - prepj->dos;
 
    if ( RandomNumber() > exp(aux) ) {
        //printf("Out swaprep\n");
        return 0; 
    }
 
    type =analyse(slice, prepj, path.nslices, path.initial_state);
    if (type==0) {
        //printf("Out swaprep\n");
        return 0;
    }

    prepj->pathlen = prepi->pathlen;
    prepj->type =  prepi->type;
    SWAP(replica[irep]->swapindex,replica[jrep]->swapindex,i);
    path.current_replica=jrep;

    path.block_stats[path.initial_state-1][prepi->index].mcacc[1].acc++;
    //printf("Out swaprep\n");

    return 1;
}




int swap_replica_0(int irep, int jrep) {

    Replica *prep;
    Slice *trial0;
    Slice *trial1;
    int i,j,pathlen0,pathlen1,inA,in1,index,maxlength,start,len,type,type0,type1,reversal;
    int final_state;
    //printf("In swaprep minus\n");

    trial0=fwtrial;
    trial1=bwtrial;

    if (irep ==0) {
        prep =replica[0];
        type =analyse_state_i(slice,prep,path.nslices, path.initial_state);  
        if (prep->type!=type)  {
            printf("swap 0->1 error : type %d and replica type %d do not match\n", type, prep->type);
            return 0;
        }
        if(type==0) {
            printf("error: replica 0 has wrong type\n");
            return 0;
        }
        reversal = (RandomNumber()<0.5) ? 1 : 0;
        
        if (reversal) {
            pathlen0 = path.nslices;
            for(j=0;j<path.nslices;j++) {  
                trial1[j]= slice[path.nslices - 1 -  j];// use trial1 for temp storage, to not affect slice
            }
        } 
        else {
            for(i=0;i< path.nslices;i++) {
                trial1[i]= slice[i];// use trial1 for temp storage, to not affect slice
            }
        }

        for(i=path.nslices-1;i>=0;i--) {
            if (in_state(&trial1[i])) {
                start = i;
                break;
            }
        }

        for(i=start;i<path.nslices;i++) {
            trial0[i-start]=trial1[i];
        }

        index = path.nslices - start - 1;
        maxlength = MAXSLICES-index;
        len = trajectory(&trial0[index],maxlength);
        if (len==0) {
            printf(" Swap replica 0: trajectory too long, irep %d jrep %d initial_state %d\n",irep,jrep,path.initial_state);
            return 0;
        }

        pathlen0 = index +len+1;         
        if (pathlen0> MAXSLICES) {
            printf("error: pathlen = %d\n",pathlen0);
            return 0;
        }
        if (pathlen0 < 3) {
            printf("Path too short out minus move: irep %d jrep %d\n",irep,jrep);
            return 0;
        }

        type= analyse(trial0,replica[1],pathlen0,path.initial_state );
        if (type==0) {
            //printf("Out swaprep minus\n");
            return 0;
        }

        for(i=0;i<pathlen0;i++) {
            slice[i]=trial0[i];
        }

        replica[1]->pathlen= pathlen0;
        replica[1]->type = type;
        path.current_replica=jrep;
        path.nslices = pathlen0;
        prep = replica[1];
        type = analyse(slice,prep,prep->pathlen, path.initial_state);  
        if (prep->type!=type)  {
            printf("swap 0->1 : type %d and replica type %d do not match\n", type, prep->type);
            return 0;
        }
        if(type==0) {
            printf("swap 0->1: replica 0 has wrong type\n");
            return 0;
        }

        final_state=in_state(&slice[path.nslices-1]);
        path.final_state = final_state;
    } 


    if (jrep==0) {
        prep =replica[1];
        type =analyse(slice, prep, path.nslices, path.initial_state);  
        if (prep->type!=type) {
            printf("error 1->0 : type %d and replica type %d do not match\n", type, prep->type);
            return 0;
        }

        if(type==0) {
            //printf("reject: replica 1 has wrong type\n");
            //printf("Out swaprep minus\n");
            return 0;
        }
    
        reversal = (RandomNumber()<0.5) ? 1 : 0;
        if ((type == 2) && (reversal ==0)) {
            reversal = 1;
            //printf("type is 2, and reversal =0, still accept\n");
            //return 0;
        }

        if (reversal) {
            pathlen1 =path.nslices;
            for(j=0;j< pathlen1;j++) { 
                trial0[j]= slice[pathlen1 - 1 -  j]; // use trial0 for temp storage
            }
        }
        else {
            for(i=0; i<path.nslices; i++) {
                trial0[i]= slice[i];// use trial0 for temp storage, to not affect slice
            }
        }
  
        for(i=path.nslices-1;i>=0;i--) {
            if (in_upper_window(&trial0[i],replica[0],path.initial_state)) {
                start = i;
                break;
            }
        }
        
        for(i=start;i<path.nslices;i++) {
            trial1[i-start]=trial0[i];
        }

        index= path.nslices - start-1;
        maxlength= MAXSLICES-index;
        len = trajectory_state_i(&trial1[index],prep,maxlength,path.initial_state);
        if (len==0) {
            printf(" Swap replica 0: trajectory too long\n");
            return 0;
        }

        pathlen1 = index +len+1;
        if (pathlen1 > MAXSLICES) {
            //printf("error: pathlen = %d\n",pathlen1);
            //printf("Out swaprep minus\n");
            return 0;
        }
        if (pathlen1 < 3) {
            printf("Path too short in of minus move: irep %d jrep %d\n",irep,jrep);
            return 0;
        }
        
        type =analyse_state_i(trial1,replica[0],pathlen1, path.initial_state);
        if ( type==0) {
            printf(" Swap replica 0: trajectory too long\n");
            return 0;
        } 
        type1=type;

        for(i=0;i<pathlen1;i++) {
            slice[i]=trial1[i];
        }
        replica[0]->pathlen= pathlen1;
        replica[0]->type = type1;
        path.current_replica=jrep;
        path.nslices=pathlen1;
        prep = replica[0];
    
        type =analyse_state_i(slice,prep,prep->pathlen, path.initial_state);  
        if (prep->type!=type)  {
            printf("swap 1->0 : type %d and replica type %d do not match\n", type, prep->type);
            return 0;
        }
        if(type==0) {
            printf("swap 1->0: replica 0 has wrong type\n");
            return 0;
        }
    }
  
    SWAP(replica[irep]->swapindex,replica[jrep]->swapindex,i);

    path.block_stats[path.initial_state-1][prep->index].mcacc[1].acc++;
    
    //printf("Out swaprep minus\n");

    return 1;
}






int swap_states(int irep,int jrep) {

    Replica *prepi,*prepj;
    int ipart,krep,islice,initial_state,final_state,pathlen,type,start,index,maxlength,len;
    double aux;

    //printf("In swapstate\n");
    
    if(path.stateswapbias==2) {
        state[path.initial_state-1].dos+=path.scalefactor;
    }
    
    prepi=replica[irep];

    if (prepi->type!=2) {
        return 0;
    }
    

    initial_state=in_state(&slice[0]);
    final_state=in_state(&slice[path.nslices-1]);

    //if(final_state==1) {
    //    return 0;
    //}

    if(final_state==initial_state) {
        printf("Warning swap states: final state is not initial state despite being type 2\n");
        return 0;
    }

    if (initial_state != path.initial_state ) {
        printf("error: initial state not correct\n");
        return 0;
    }

    if(final_state==0) {
        printf("error state swap: final state unknown\n");
        return 0;
    }

    if(initial_state==0) {
        printf("error state swap: initial state unknown\n");
        return 0;
    }

    path.block_stats[path.initial_state-1][prepi->index].mcacc[2].tries++;

    //pick random replica of the final state;
    jrep = (int)(RandomNumber()*(state[final_state-1].nrep-1)) +1;
    if(jrep==0) {
        printf("Warning swap states: jrep is 0, should at least be 1\n");
    }
    prepj=&state[final_state-1].srep[jrep];  
    
    //jrep = state[final_state-1].nrep-1;
    //prepj=&state[final_state-1].srep[jrep];  
    
    double nrepratio;
    nrepratio = (double)(state[final_state-1].nrep-1.)/(double)(state[initial_state-1].nrep-1.);
    //nrepratio=1;
    if(path.stateswapbias==1) {
        aux = prepi->dos - prepj->dos;
        if (RandomNumber() > (nrepratio*exp(aux))) {
            return 0; 
        }
    }
    else if(path.stateswapbias==2) {
        aux = state[initial_state-1].dos-state[final_state-1].dos;
        if (RandomNumber() > (nrepratio*exp(aux))) {
            //printf("Out swapstate\n");
            return 0; 
        }
    }
    else {
        if (RandomNumber() > nrepratio) {
            //printf("Out swapstate\n");
            return 0; 
        }
    }
    //printf("Rejected state swap I %d J %d: nrepratio: %lf dosdiff: %lf\n",initial_state,final_state,nrepratio,aux);

    path.initial_state = final_state;
    path.final_state = initial_state;
    path.current_gsen = state[final_state-1].target.energy;
    for(ipart=0; ipart<sys.npart; ipart++) {
        if(slice[path.nslices-1].clusterid[ipart] == slice[path.nslices-1].indexmaxcluster) {
            path.current_maxclustid[ipart]=1;
        }
        else {
            path.current_maxclustid[ipart]=0;
        }
    }
    //path.current_gsen = state[final_state-1].target.enmaxclust;


    for(islice=0;islice< path.nslices;islice++) {
        trial[islice]= slice[path.nslices - 1 -  islice];
        create_all_rc(&(trial[islice]));
    }


    type = analyse(trial,prepj,path.nslices,final_state);
    if ( type!=2 ) {
        //printf("Rejected state swap:wrong type %d \n",type);
        //printf("Reverse path from state %d to state %d, initial_state=%d\n",initial_state,final_state,path.initial_state);
        //printf("Swap state is attempted from irep %d to jrep %d\n",irep,jrep);
        //printf("New current gsen %lf\n",path.current_gsen);
        //printf("Necessary lambda value to reach: %lf\n",state[final_state-1].lambda[jrep]);
        //printf("Maximum lambda value reached after swap: %lf\n",trial_maxop);

        //resetting path parameters
        path.current_gsen = state[initial_state-1].target.energy;
        for(ipart=0; ipart<sys.npart; ipart++) {
            if(slice[0].clusterid[ipart] == slice[0].indexmaxcluster) {
                path.current_maxclustid[ipart]=1;
            }
            else {
                path.current_maxclustid[ipart]=0;
            }
        }
        //path.current_gsen = state[initial_state-1].target.enmaxclust;
        path.initial_state = initial_state;
        path.final_state = final_state;
        //printf("Out swapstate\n");
        return 0;
    }
  
    path.nreplica = state[final_state-1].nrep;

    for(krep=0; krep<path.nreplica; krep++) {
        replica[krep] = &state[path.initial_state-1].srep[krep];
    }

    prepj->pathlen = path.nslices;
    for(islice=0; islice<path.nslices; islice++) {
        slice[islice] = trial[islice];  
    }

    type = analyse(slice,prepj, prepj->pathlen, path.initial_state);
    prepj->type = type;
    path.current_replica=jrep;
 
    path.block_stats[initial_state-1][prepi->index].mcacc[2].acc++;

    //printf("Out swapstate\n");

    return 1;
}



int reverse_replica(int irep) {

    Replica *prep;
    int islice,type,ipart,maxclustchanged;

    //printf("In repreverse\n");

    prep=replica[irep];
    path.block_stats[path.initial_state-1][irep].mcacc[3].tries++;

    if (irep==0) {
        //printf("Out repreverse\n");
        return 0;  
    }

    type = analyse(slice,prep, path.nslices, path.initial_state);
    if (prep->type != type) {
        printf("error swap rep: type %d and replica type %d do not match\n",type,prep->type);
    }

    if (prep->type!=1) {
        //printf("Out repreverse\n");
        return 0; //reversal can only happen for trajectories who end and start in the same state
    }


    maxclustchanged=0;
    for(ipart=0; ipart<sys.npart; ipart++) {
        if(slice[path.nslices-1].clusterid[ipart] == slice[path.nslices-1].indexmaxcluster) {
            if(path.current_maxclustid[ipart]!=1) {
                maxclustchanged=1;
            }
            else {
                maxclustchanged=0;
            }
            path.current_maxclustid[ipart]=1;
        }
        else {
            path.current_maxclustid[ipart]=0;
        }
    }

    //for(ipart=0; ipart<sys.npart; ipart++) {
    //    if(slice[path.nslices-1].clusterid[ipart] == slice[path.nslices-1].indexmaxcluster) {
    //        path.current_maxclustid[ipart]=1;
    //    }
    //    else {
    //        path.current_maxclustid[ipart]=0;
    //    }
    //}


    for(islice=0;islice< prep->pathlen;islice++) {
        trial[islice]= slice[prep->pathlen - 1 -  islice];
        if(maxclustchanged==1) {
            create_all_rc(&trial[islice]);
        }
    }



    type = analyse(trial,prep,prep->pathlen,path.initial_state);
    if ( type==0) {
        return 0;
    }
    if ( type==2) {
        return 0;
    }

    for(islice=0;islice< prep->pathlen;islice++) {
        slice[islice] = trial[islice];
    }

  
    path.block_stats[path.initial_state-1][irep].mcacc[3].acc++;

    //printf("Out repreverse\n");

    return 1;
}



int add_newstate(Slice *psl) {

    vector rcom,dr;
    double r2,nexteps,maxeps;
    int ipart,jpart,istate,j,k,islice;

    if(path.nstates>=MAXSTATES) {
        printf("Warning: too many states\n");
        printf("Did not add possible new state\n");
        return 0;
    }


    if(psl->bonddidnotbreak==1) {
        printf("Warning: not all bonds completely broken... Can not be a new state.\n");
        return 0;
    }

    if(psl->nbonds==0) {
        printf("Warning: no bonds formed... Will lead to unbound state.\n");
        return 0;
    }


    istate=path.nstates;
    psl->energy = potential_energy(psl);
    print_rc(psl,path.nstates);
    state[istate].target=*psl;

    //setting replicas of new state based on energy
    state[istate].volume_op = -0.3*state[istate].target.energy;
    //state[istate].volume_op = -0.3*state[istate].target.enmaxclust;
    //state[istate].volume_opdist=1000000000;
    state[istate].srep[0].lambda=state[istate].volume_op;
    state[istate].srep[1].lambda=state[istate].volume_op;
    printf("For state %d rep 0 has lambda %lf\n",istate,state[istate].srep[0].lambda);
    printf("For state %d rep 1 has lambda %lf\n",istate,state[istate].srep[1].lambda);
    state[istate].nrep=2;
    //maxeps = -state[istate].target.enmaxclust;
    maxeps = -state[istate].target.energy;
    nexteps = state[istate].srep[1].lambda;
    do {
        nexteps+=path.diffenrep;
        if(nexteps>(maxeps-0.1)) {
            break;
        }
        else if(state[istate].nrep==MAXREPLICA) {
            break;
        }
        else {
            state[istate].srep[state[istate].nrep].lambda=nexteps;
            printf("For state %d rep %d has lambda %lf\n",istate,state[istate].nrep,state[istate].srep[state[istate].nrep].lambda);
            state[istate].nrep++;
        }
    } while(TRUE);
    //printf("Total replicas defined for state %d is %d\n",istate,state[istate].nrep);

    path.nstates++;
    if(check_stability(istate)==0) {
        printf("Proposed state not stable enough: rejection state\n");
        path.nstates--;
        return 0;
    }

    printf("Added state %d to database. Now %d number of states defined.\n",istate,path.nstates);

    return 1;
}


int check_stability(int istate) {

    int knownstate,islice,trial_state,relaxed2other,i;
    double oldmobr,oldmobrnoise,avsumopdist,countinstate;

    trial[0] = state[istate].target;
    trial_state=istate+1;
    relaxed2other=0;
    //countinstate=0;
    //avsumopdist=0;
    for(islice=1; islice<MAXSLICES; islice++) {
        trial[islice]=trial[islice-1];
        propagate_bd(&trial[islice]);
        create_all_rc(&trial[islice]);
        i=in_state(&trial[islice]);
        if(i!=trial_state && i!=0) {
            knownstate=i;
            relaxed2other=1;
            break;
        }
        //if(i==trial_state) {
        //    avsumopdist+=trial[islice].sumopdist[istate];
        //    countinstate++;
        //}
    }

    if(relaxed2other) {
        printf("State is not stable. It has relaxed to state %d after %d slices.\n",knownstate,islice);
        return 0;
    }

    //avsumopdist/=countinstate;
    //state[istate].volume_opdist=avsumopdist;
    
    printf("State is found to be stable enough\n");
    return 1;
}












