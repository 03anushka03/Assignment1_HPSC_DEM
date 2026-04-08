#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <chrono>

using namespace std;

const double DT=1.0e-4;
const double TOTAL_TIME=2.0;
const double KN=1.0e5;
const double GAMMA_N=10.0;
const double GX=0.0;
const double GY=0.0;
const double GZ=-9.81;
const double LX=1.0;
const double LY=1.0;
const double LZ=1.0;
const double RADIUS=0.02;
const double MASS=1.0;
const int NMAX=6000;

double px[NMAX], py[NMAX], pz[NMAX];
double vx[NMAX], vy[NMAX], vz[NMAX];
double fx[NMAX], fy[NMAX], fz[NMAX];
double radius[NMAX], mass_arr[NMAX];

int N=0;


void initialise(int n, double x0, double y0, double z0, double u0, double v0, double w0, bool random_pos);
void zero_forces();
void add_gravity();
void particle_contacts();
void wall_contacts();
void integrate_step();
double kinetic_energy();
void write_snapshot(const string& fname);

int main(int argc, char* argv[])
{
  string mode="multi";
  int npart=200;
  if (argc>1)
  {
     mode=argv[1];
  }
  if (argc>2){
     npart=atoi(argv[2]);

   double sim_time=TOTAL_TIME;
   bool do_pp_contacts=true;
   bool do_walls=true;
   bool do_gravity=true;

   if(mode=="freefall")
   {
      N=1;
      initialise(1,0.5,0.5,0.9,0.0,0.0,0.0, false);
      do_pp_contacts=false;
      sim_time=1.0;
      cout <<"MODE: Free Fall\n"
           <<" ball starts at z=0.9 with zero velocity\n"
           <<" Analytical: z(t)=0.9-4.905*t^2\n\n";
    }
   else if(mode=="constvel")
   {
       N=1;
       initialise(1, 0.1, 0.5, 0.5, 2.0, 1.0, 0.5, false);
       do_pp_contacts=false;
       do_walls=false;
       do_gravity=false;
       sim_time=0.4;
       cout<<"MODE: Constant Velocity (gravity OFF)\n"
           <<" Initial velocity: vx2.0, vy=1.0, vz=0.5m/s\n";
   }
   else if(mode=="bounce")
   {
       N=1;
       initialise(1,0.5,0.5,0.3,0.0,0.0,0.0, false);
       do_pp_contacts=false;
       sim_time=2.0;
       cout<<"MODE: Bouncing Ball\n"
           <<" Ball starts at z=0.3, drops and bounces\n"
           <<" Expect: each rebound is lower than the last\n\n";
   }
   else
   {
       N=npart;
       initialise(N, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, true);
       cout <<"MODE: Multi-particle simulation\n"
            <<" N= "<< N << " particles\n\n"; 
   }

   cout << "Total steps: "<< (int)(sim_time/DT)
        <<" | dt= "<< DT <<" s\n\n";
   
   ofstream f_traj(mode+ "_trajectory.dat");
   ofstream f_ke (mode+"_energy.dat");
   f_traj << "# time x y z vz\n";
   f_ke << "# time kinetic_energy\n";
   
   if(mode=="freefall")
   {
     ofstream f_anal("freefall_analytical.dat");
     f_anal << "#time z_analytical vz_analytical\n";
     for(double t=0.0; t<=sim_time; t+=DT)
     {
        double z_a=0.9+0.5*GZ*t*t;
        double vz_a=GZ*t;
        f_anal << t << " " << z_a << " " << vz_a << "\n";
     }
   }

   auto t_start = chrono::high_resolution_clock::now();
 
    double tm_gravity  = 0.0;
    double tm_contacts = 0.0;
    double tm_walls    = 0.0;
    double tm_integ    = 0.0;
 
    auto now     = []{ return chrono::high_resolution_clock::now(); };
    auto elapsed = [](auto a, auto b){
        return chrono::duration<double>(b - a).count();
    };

    double time=0.0;
    int step=0;

    while(time<sim_time)
    {
      zero_forces();

      auto t0=now();
      if(do_gravity) add_gravity();
      tm_gravity+=elapsed(t0, now());

      t0=now();
      if(do_pp_contacts) particle_contacts();
      tm_contacts+=elapsed(t0, now());

      t0=now();
      if(do_walls) wall_contacts();
      tm_walls+=elapsed(t0, now());

      t0=now();
      integrate_step();
      tm_integ+=elapsed(t0, now());

      if(step%100==0)
      {
         double KE=kinetic_energy();
         f_ke << time << " " << KE << "\n";
         f_traj << time << " " << px[0]<< " " << py[0] << pz[0] << " "<< vz[0]<< "\n";      
      }

      if(step%10000==0)
        write_snapshot(mode+"_snap_"+to_string(step)+".dat");

      time+=DT;
      step++;
      }

      double t_total=elapsed(t_start, now());

   cout << "=== Simulation complete ===\n";
   cout << "Total wall-clock time : " << t_total << "s\n\n";
   cout << "--- Time per function ---\n";
   cout << "Gravity      : "<< tm_gravity << "s (" << 100.0*tm_gravity/t_total << " %)\n";
   cout << " P-P contacts L "<< tm_contacts << "s (" << 100.0*tm_contacts/t_total << "%)\n";
   cout << " Wall contacts : "<< tm_walls << "s (" << 100.0*tm_walls/t_total << "%)\n";
   cout << " Integration : "<< tm_integ << "s ("<< 100.0*tm_integ/t_total<< " %)\n\n";
      
   ofstream f_time("timing_serial_"+mode+ ".dat", ios::app);
   f_time << N << " " << t_total << "\n";

   return 0; 
}
}


void initialise(int n, double x0, double y0, double z0, double u0, double v0, double w0, bool random_pos)
{
   srand48(42);
   for(int i=0;i<n;i++)
   {
      radius[i]=RADIUS;
      mass_arr[i]=MASS;

      if(random_pos)
      {
         double margin=RADIUS;
         px[i]=margin+drand48()+(LX-2.0*margin);
         py[i]=margin+drand48()*(LY-2.0*margin);
         pz[i]=margin+drand48()*(LZ-2.0*margin);
      }
      else
      {
         px[i]=x0;
         py[i]=y0;
         pz[i]=z0;
      }

      vx[i]=u0;
      vy[i]=v0;
      vz[i]=w0;

      fx[i]=0.0;
      fy[i]=0.0;
      fz[i]=0.0;
   }
}


void zero_forces()
{
   for(int i=0; i<N; i++)
   {
      fx[i]=fy[i]=fz[i]=0.0;
   }
}

void add_gravity()
{
   for(int i=0; i<N; i++)
   {
      fx[i]+=mass_arr[i]*GX;
      fy[i]+=mass_arr[i]*GY;
      fz[i]+=mass_arr[i]*GZ;
   }
}

void particle_contacts()
{
   for(int i=0; i<N; i++)
   {
      for(int j=i+1;j<N;j++)
      {
         double rx=px[j]=px[i];
         double ry=py[j]=py[i];
         double rz=pz[j]-pz[i];

         double dij=sqrt(rx*rx+ry*ry+rz*rz);

         double delta=radius[i]+radius[j]-dij;

         if(delta<=0.0)continue;

         double nx=rx/dij;
         double ny=ry/dij;
         double nz=rz/dij;

         double vn=(vx[j]-vx[i])*nx + (vy[j]-vy[i])*ny + (vz[j]-vz[i])*nz;

         double Fn=KN*delta-GAMMA_N*(vn>0.0?vn:0.0);
         if(Fn<0.0) Fn=0.0;

         double Fx=Fn*nx;
         double Fy=Fn*ny;
         double Fz=Fn*nz;

         fx[i]-=Fx; fy[i]-=Fy; fz[i]-=Fz;
         fx[j]+=Fx; fy[j]+=Fy; fz[j]+=Fz;
      }
   }
}

void wall_contacts()
{
   for(int i=0; i<N; i++)
   {
      double delta, vn, Fn;
      delta=radius[i]=pz[i];
      if(delta>0.0)
      {
         vn=-vz[i];
         Fn=KN*delta-GAMMA_N*(vn>0.0? vn:0.0);
         if(Fn<0.0)Fn=0.0;
         fz[i]+=Fn;
      }

      delta=(pz[i]+radius[i])-LZ;
      if(delta>0.0)
      {
         vn=vz[i];
         Fn=KN*delta-GAMMA_N*(vn>0.0?vn:0.0);
         if(Fn<0.0)Fn=0.0;
         fz[i]-=Fn;
      }

      delta=(px[i]+radius[i])-LX;
      if(delta>0.0)
      {
         vn=vx[i];
         Fn=KN*delta-GAMMA_N*(vn>0.0?vn:0.0);
         if(Fn<0.0)Fn=0.0;
         fx[i]-=Fn;
      }

      delta=radius[i]-px[i];
      if(delta>0.0)
      {
         vn=-vx[i];
         Fn=KN*delta-GAMMA_N*(vn>0.0?vn:0.0);
         if(Fn<0.0)Fn=0.0;
         fx[i]+=Fn;
      }

      delta=(py[i]+radius[i])-LY;
      if(delta>0.0)
      {
         vn=vy[i];
         Fn=KN*delta-GAMMA_N*(vn>0.0?vn:0.0);
         if(Fn<0.0)Fn=0.0;
         fy[i]-=Fn;
      }

      delta=radius[i]-py[i];
      if(delta>0.0)
      {
         vn=-vy[i];
         Fn=KN*delta-GAMMA_N*(vn>0.0?vn:0.0);
         if(Fn<0.0)Fn=0.0;
         fy[i]+=Fn;
      }
   }
}

void integrate_step()
{
   for(int i=0;i<N;i++)
   {
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
   double KE=0.0;
   for(int i=0;i<N;i++)
   {
      double v2=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
      KE+=0.5*mass_arr[i]*v2;
   }
   return KE;
}
 
void write_snapshot(const string& fname)
{
   ofstream f(fname);
   f << "# x y z vx vy vz\n";
   for(int i=0;i<N;i++)
   {
      f << px[i] << " " << py[i] << " " << pz[i] << " " << "  " << vx[i] << " " << vy[i] << " " << vz[i] << "\n";
   }
}
