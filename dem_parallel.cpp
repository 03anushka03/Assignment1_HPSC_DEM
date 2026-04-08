#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <chrono>
#include <omp.h>

using namespace std;

const double DT=1e-4,TOTAL_TIME=2.0;
const double KN=1e5,GAMMA_N=10.0;
const double GX=0.0,GY=0.0,GZ=-9.81;
const double LX=1.0,LY=1.0,LZ=1.0;
const double RADIUS=0.02,MASS=1.0;
const int NMAX=6000;

double px[NMAX],py[NMAX],pz[NMAX];
double vx[NMAX],vy[NMAX],vz[NMAX];
double fx[NMAX],fy[NMAX],fz[NMAX];
double radius[NMAX],mass_arr[NMAX];

int N=0;

void initialise(int,double,double,double,double,double,double,bool);
void zero_forces();
void add_gravity();
void particle_contacts();
void wall_contacts();
void integrate_step();
double kinetic_energy();
void write_snapshot(const string&);

int main(int argc,char* argv[])
{
    string mode=(argc>1)?argv[1]:"multi";
    int npart=(argc>2)?atoi(argv[2]):200;

    cout<<"Threads: "<<omp_get_max_threads()<<"\n";

    double sim_time=TOTAL_TIME;
    bool do_pp=true,do_walls=true,do_gravity=true;

    if(mode=="freefall"){
        N=1;
        initialise(1,0.5,0.5,0.9,0,0,0,false);
        do_pp=false;do_walls=false;
        sim_time=1.0;
    }
    else if(mode=="constvel"){
        N=1;
        initialise(1,0.1,0.5,0.5,2,1,0.5,false);
        do_pp=false;do_walls=false;do_gravity=false;
        sim_time=0.4;
    }
    else if(mode=="bounce"){
        N=1;
        initialise(1,0.5,0.5,0.3,0,0,0,false);
        do_pp=false;
    }
    else{
        N=npart;
        initialise(N,0,0,0,0,0,0,true);
    }

    ofstream f_traj(mode+"_trajectory.dat");
    ofstream f_ke(mode+"_energy.dat");

    auto t_start=chrono::high_resolution_clock::now();
    auto now=[](){return chrono::high_resolution_clock::now();};
    auto elapsed=[](auto a,auto b){return chrono::duration<double>(b-a).count();};

    double time=0;
    int step=0;

    while(time<sim_time){
        zero_forces();

        if(do_gravity) add_gravity();
        if(do_pp) particle_contacts();
        if(do_walls) wall_contacts();

        integrate_step();

        if(step%100==0){
            f_ke<<time<<" "<<kinetic_energy()<<"\n";
            f_traj<<time<<" "<<px[0]<<" "<<py[0]<<" "<<pz[0]<<" "<<vz[0]<<"\n";
        }

        if(step%10000==0)
            write_snapshot(mode+"_snap_"+to_string(step)+".dat");

        time+=DT;
        step++;
    }

    cout<<"Total time: "<<elapsed(t_start,now())<<" s\n";
    return 0;
}

void initialise(int n,double x0,double y0,double z0,double u0,double v0,double w0,bool random_pos)
{
    srand48(42);
    for(int i=0;i<n;i++){
        radius[i]=RADIUS;
        mass_arr[i]=MASS;

        if(random_pos){
            double m=RADIUS;
            px[i]=m+drand48()*(LX-2*m);
            py[i]=m+drand48()*(LY-2*m);
            pz[i]=m+drand48()*(LZ-2*m);
        }else{
            px[i]=x0;py[i]=y0;pz[i]=z0;
        }

        vx[i]=u0;vy[i]=v0;vz[i]=w0;
        fx[i]=fy[i]=fz[i]=0.0;
    }
}

void zero_forces()
{
#pragma omp parallel for
    for(int i=0;i<N;i++) fx[i]=fy[i]=fz[i]=0.0;
}

void add_gravity()
{
#pragma omp parallel for
    for(int i=0;i<N;i++){
        fx[i]+=mass_arr[i]*GX;
        fy[i]+=mass_arr[i]*GY;
        fz[i]+=mass_arr[i]*GZ;
    }
}

void particle_contacts()
{
#pragma omp parallel for schedule(dynamic,16)
    for(int i=0;i<N;i++){
        for(int j=i+1;j<N;j++){
            double rx=px[j]-px[i];
            double ry=py[j]-py[i];
            double rz=pz[j]-pz[i];

            double dij=sqrt(rx*rx+ry*ry+rz*rz);
            double delta=radius[i]+radius[j]-dij;
            if(delta<=0) continue;

            double nx=rx/dij,ny=ry/dij,nz=rz/dij;
            double vn=(vx[j]-vx[i])*nx+(vy[j]-vy[i])*ny+(vz[j]-vz[i])*nz;

            double Fn=KN*delta-GAMMA_N*(vn>0?vn:0);
            if(Fn<0) Fn=0;

            double Fx=Fn*nx,Fy=Fn*ny,Fz=Fn*nz;

#pragma omp atomic
            fx[i]-=Fx;
#pragma omp atomic
            fy[i]-=Fy;
#pragma omp atomic
            fz[i]-=Fz;

#pragma omp atomic
            fx[j]+=Fx;
#pragma omp atomic
            fy[j]+=Fy;
#pragma omp atomic
            fz[j]+=Fz;
        }
    }
}

void wall_contacts()
{
#pragma omp parallel for
    for(int i=0;i<N;i++){
        double delta,vn,Fn;

        delta=radius[i]-pz[i];
        if(delta>0){vn=-vz[i];Fn=KN*delta-GAMMA_N*(vn>0?vn:0);if(Fn<0)Fn=0;fz[i]+=Fn;}

        delta=(pz[i]+radius[i])-LZ;
        if(delta>0){vn=vz[i];Fn=KN*delta-GAMMA_N*(vn>0?vn:0);if(Fn<0)Fn=0;fz[i]-=Fn;}

        delta=radius[i]-px[i];
        if(delta>0){vn=-vx[i];Fn=KN*delta-GAMMA_N*(vn>0?vn:0);if(Fn<0)Fn=0;fx[i]+=Fn;}

        delta=(px[i]+radius[i])-LX;
        if(delta>0){vn=vx[i];Fn=KN*delta-GAMMA_N*(vn>0?vn:0);if(Fn<0)Fn=0;fx[i]-=Fn;}

        delta=radius[i]-py[i];
        if(delta>0){vn=-vy[i];Fn=KN*delta-GAMMA_N*(vn>0?vn:0);if(Fn<0)Fn=0;fy[i]+=Fn;}

        delta=(py[i]+radius[i])-LY;
        if(delta>0){vn=vy[i];Fn=KN*delta-GAMMA_N*(vn>0?vn:0);if(Fn<0)Fn=0;fy[i]-=Fn;}
    }
}

void integrate_step()
{
#pragma omp parallel for
    for(int i=0;i<N;i++){
        vx[i]+=(fx[i]/mass_arr[i])*DT;
        vy[i]+=(fy[i]/mass_arr[i])*DT;
        vz[i]+=(fz[i]/mass_arr[i])*DT;

        px[i]+=vx[i]*DT;
        py[i]+=vy[i]*DT;
        pz[i]+=vz[i]*DT;
    }
}

double kinetic_energy()
{
    double KE=0;
    for(int i=0;i<N;i++){
        double v2=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
        KE+=0.5*mass_arr[i]*v2;
    }
    return KE;
}

void write_snapshot(const string& fname)
{
    ofstream f(fname);
    for(int i=0;i<N;i++)
        f<<px[i]<<" "<<py[i]<<" "<<pz[i]<<" "
         <<vx[i]<<" "<<vy[i]<<" "<<vz[i]<<"\n";
}