#include <complex.h>
#include <fftw3.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <stdlib.h>
#include "mpi.h"

//#include "angf.h"

using namespace std;


int main(int argc, char* argv[]) {
  
  int mpisize, mpirank;

  MPI::Init(argc,argv);
  mpisize=MPI::COMM_WORLD.Get_size();
  mpirank=MPI::COMM_WORLD.Get_rank();


  vector<double> tt; // timing 
  double tm=0.0;
  double DT=0.0;
  double DT_max=0.0;

  MPI::COMM_WORLD.Barrier();
  tm=MPI::Wtime();
  tt.push_back(tm);
  MPI::COMM_WORLD.Barrier();

  int ndims=2;
  int dims[2];
  dims[0]=atoi(argv[1]); // proc # in row: receiver index
  dims[1]=atoi(argv[2]); // proc # in col: time segment index
  if (dims[0]*dims[1]!=mpisize) {
    cerr << "dims[0]*dims[1] not equal to mpisize!!! exiting" << endl;
    exit(1);
  }

  bool periods[2]={true,false};
  bool reorder=true;

  MPI::Cartcomm com2d, comcol, comrow;
  
  com2d=MPI::COMM_WORLD.Create_cart(ndims,dims,periods,reorder);
  int coord[2];
  com2d.Get_coords(mpirank,ndims,coord);


  bool remaincol[2]={true,false};
  bool remainrow[2]={false,true};

  comcol=com2d.Sub(remaincol);
  comrow=com2d.Sub(remainrow);

  int colrank=comcol.Get_rank();
  int rowrank=comrow.Get_rank();
  /*
  ofstream fou("psin.ou"+to_string(mpirank));
  fou << mpirank << " " << colrank << " " << rowrank << " " << coord[0] << " " << coord[1] << endl;
  fou.close();
  */



  //cout << "rank: " << mpirank << endl;

  MPI::COMM_WORLD.Barrier();
  tm=MPI::Wtime();
  DT=tm-tt.back();
  MPI::COMM_WORLD.Reduce(&DT,&DT_max,1,MPI_DOUBLE,MPI_MAX,0);
  if(mpirank==0) cout << "setup cart elapse time: " << DT_max  << " sec" << endl;
  tt.push_back(tm);
  MPI::COMM_WORLD.Barrier();


  ////


  int sps=atoi(argv[4]);
  double dt=1.0/sps;

  int Nr=0;
  int Nr_row=0;
  vector<string> binlst_row;
  vector<int> rlst_row;
 
  if (rowrank==0) {
    ifstream finp(argv[3]);
    if(!finp.is_open()){
      cerr << "error opening file: " << argv[3] << endl;
      exit(1);
    }
    vector<int> rlst;
    vector<string> binlst;
    string line;
    int ir=0;
    while(getline(finp,line)) {
      binlst.push_back(line);
      rlst.push_back(ir);
      ++ir;
    }
    finp.close();
    //cout << "rank: " << mpirank << " " << rlst.size() << endl;
    
    Nr=binlst.size();
    
    vector<int> NrPerRow(dims[0],0);
    ir=0;
    for(int i=0; i!=Nr; ++i) {
      NrPerRow[ir] += 1;
      if(ir==dims[0]-1)
	ir=0;
      else
	++ir;
    }
    //cout << mpirank << " " << NrPerRow[0] << " " << NrPerRow[1] << endl;
    Nr_row=NrPerRow[colrank];

    vector<int> binB(dims[0],0);
    for(int i=1; i!=dims[0]; ++i)
      binB[i]=binB[i-1]+NrPerRow[i-1];

    for(int i=0; i!=Nr_row; ++i) {
      binlst_row.push_back(binlst[binB[colrank]+i]);
      rlst_row.push_back(rlst[binB[colrank]+i]);
    }

    /*
    if(mpirank==8) {
      for(int i=0; i!=Nr_row; ++i)
	cout << binlst_row[i] << endl;
    }
    */

  }

  comrow.Bcast(&Nr, 1, MPI_INT, 0);  
  comrow.Bcast(&Nr_row, 1, MPI_INT, 0);  

  rlst_row.resize(Nr_row);
  comrow.Bcast(&rlst_row[0], Nr_row, MPI_INT, 0);

  //cout << mpirank << " " << Nr << " " << Nr_row << " " << rlst_row.size() <<  endl;
  /*
  if(colrank==15) {
    cout << "rank: " << mpirank << " ";
    for(int i=0; i!=rlst_row.size(); ++i)
      cout << rlst_row[i] << " ";
    cout <<endl;
  }
  
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Abort(1);
  */

  MPI::COMM_WORLD.Barrier();
  tm=MPI::Wtime();
  DT=tm-tt.back();
  MPI::COMM_WORLD.Reduce(&DT,&DT_max,1,MPI_DOUBLE,MPI_MAX,0);
  if(mpirank==0) cout << "proc psin.in elapse time: " << DT_max << " sec" << endl;
  tt.push_back(tm);
  MPI::COMM_WORLD.Barrier();


  ////


  long NN=0;
  double astime=0, aetime=0;

  if(rowrank==0) {

    double astimel=1e30;
    double aetimel=-1e30;

    for (int ir=0; ir!=Nr_row; ++ir) {

      ifstream binp(binlst_row[ir].c_str(),ios::in|ios::binary|ios::ate);
      if(!binp.is_open()){
	cerr << "error opening file: " << binlst_row[ir] << endl;
	exit(1);
      }
      
      size_t N  = (size_t)binp.tellg()-3*sizeof(double)-2*sizeof(float);
      size_t NS = N/sizeof(float);
      binp.seekg(0,binp.beg);
      double t0=0.0;
      binp.read(reinterpret_cast<char*>(&t0),sizeof(double));
      binp.close();
      double t1=t0+(NS-1)*dt;

      if(t0<astimel) astimel=t0;
      if(t1>aetimel) aetimel=t1;
    }

    comcol.Allreduce(&astimel,&astime,1,MPI_DOUBLE,MPI_MIN);
    comcol.Allreduce(&aetimel,&aetime,1,MPI_DOUBLE,MPI_MAX);

    NN=round((aetime-astime)/dt)+1;
    cout << mpirank << scientific << setprecision(10) << " astime = " << astime << " aetime = " << aetime << " " << Nr_row << " " << NN << endl;
    
  } // if rowrank==0

  /*
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Abort(1);
  */

  comrow.Bcast(&NN,1,MPI_LONG,0);
  long Ns=atol(argv[5]);
  int  NsegPerCol=floor((double)NN/(double)Ns/(double)dims[1]);
  long NsPerCol=NsegPerCol*Ns;
  long NN1=NsPerCol*dims[1];

  cout << "rank: " << mpirank << " " << Nr_row << " " << Ns << " " << NsegPerCol << " " << NsPerCol << " " << NN1 << " " << NN << endl;


  MPI::COMM_WORLD.Barrier();
  tm=MPI::Wtime();
  DT=tm-tt.back();
  MPI::COMM_WORLD.Reduce(&DT,&DT_max,1,MPI_DOUBLE,MPI_MAX,0);
  if(mpirank==0) cout << "get array time elapse time: " << DT_max << " sec" << endl;
  tt.push_back(tm);
  MPI::COMM_WORLD.Barrier();


  ////


  float *S;
  float *ss=new float[Nr_row*NsPerCol];

  if(rowrank==0)
    S=new float[NN1]();
  
  double tsc0=0.0;
  double tsc1=0.0;
  double tscs=0.0;

  for (int ir=0; ir!=Nr_row; ++ir) {

    if(rowrank==0) {

      ifstream binp(binlst_row[ir].c_str(),ios::in|ios::binary|ios::ate);
      size_t N=(size_t)binp.tellg()-3*sizeof(double)-2*sizeof(float);
      size_t NS = N/sizeof(float);
      binp.seekg(0,binp.beg);
    
      double t0, xc, yc;
      float aID, rID;
      float* s=new float[NS];
      
      binp.read(reinterpret_cast<char*>(&t0),sizeof(double));
      binp.read(reinterpret_cast<char*>(&xc),sizeof(double));
      binp.read(reinterpret_cast<char*>(&yc),sizeof(double));
      binp.read(reinterpret_cast<char*>(&aID),sizeof(float));
      binp.read(reinterpret_cast<char*>(&rID),sizeof(float));
      binp.read(reinterpret_cast<char*>(s), N);
      binp.close();
    
      float* z=new float[NN]();
      //for(int i=0; i!=NN; ++i) z[i]=0.0;
      size_t i0=round((t0-astime)/dt);
      for(size_t i=0; i!=NS; ++i) z[i+i0]=s[i];
      delete[] s;
      
      for(size_t i=0; i!=NN1; ++i) S[i]=z[i];
      delete[] z;
      /*
      if(colrank==16) {
	cout << binlst_row[ir] << " " << i0 << endl;
      }
      */

    } // if rowrank==0

    //MPI::COMM_WORLD.Barrier(); // wrong logic if different rows have different number of receivers
    comrow.Barrier();
    tsc0=MPI::Wtime();
    comrow.Scatter(S,NsPerCol,MPI_FLOAT,&ss[ir*NsPerCol],NsPerCol,MPI_FLOAT,0);
    tsc1=MPI::Wtime();
    tscs += (tsc1-tsc0);
    if(rowrank==0) {
      cout << "colrank: " << colrank << " ir = " << ir << " tscs = " << tscs << endl;
    }
    //MPI::COMM_WORLD.Barrier(); // wrong logic if different rows have different number of receivers
    comrow.Barrier();
    
  } // for ir

  if(rowrank==0)
    delete[] S;

  /*
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Abort(1);
  */

  /*
  if (colrank==12) {
    stringstream iss;
    iss << setw(2) << setfill('0') << rowrank;
    string colid; iss >> colid;
    ofstream fou(("psin.ou00"+colid).c_str(), ios::out|ios::binary);
    fou.write(reinterpret_cast<char*>(&ss[NsPerCol]),sizeof(float)*(NsPerCol));
    fou.close();
  }
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Abort(1);
  */

  MPI::COMM_WORLD.Barrier();
  tm=MPI::Wtime();
  DT=tm-tt.back();
  MPI::COMM_WORLD.Reduce(&DT,&DT_max,1,MPI_DOUBLE,MPI_MAX,0);
  double tscs_max=0.0;
  MPI::COMM_WORLD.Reduce(&tscs,&tscs_max,1,MPI_DOUBLE,MPI_MAX,0);
  if(mpirank==0) {
    cout << "read/scatter recordings elapse time: " << DT_max << " sec" << endl;
    cout << "scatter elapse time: " << tscs_max << " sec" << endl;
  }
  tt.push_back(tm);
  MPI::COMM_WORLD.Barrier();

  
  //// time domain normalization


  float wins=atof(argv[6]);
  int winp=round(wins*sps);
  //cout <<"rank: " << mpirank << " " << winp << endl;

  double* tc=new double[NsPerCol];
  double* wt=new double[NsPerCol];

  for(int ir=0; ir!=Nr_row; ++ir) {

    tc[0]=fabs(ss[ir*NsPerCol]);
    for(size_t i=1; i!=NsPerCol; ++i)
      tc[i]=tc[i-1]+fabs(ss[ir*NsPerCol+i]);
    
    for(size_t i=winp; i!=NsPerCol-winp; ++i){
      double tmp=tc[i+winp]-tc[i-winp];
      wt[i]=(tmp==0?1:tmp);
    }
    
    for(size_t i=winp; i!=NsPerCol-winp; ++i)
      ss[ir*NsPerCol+i]/=wt[i];
  }

  delete[] tc;
  delete[] wt;

  /*
  if (colrank==12) {
    stringstream iss;
    iss << setw(2) << setfill('0') << rowrank;
    string colid; iss >> colid;
    ofstream fou(("psin.ou00"+colid).c_str(), ios::out|ios::binary);
    fou.write(reinterpret_cast<char*>(&ss[NsPerCol]),sizeof(float)*(NsPerCol));
    fou.close();
  }
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Abort(1); 
  */

  MPI::COMM_WORLD.Barrier();
  tm=MPI::Wtime();
  DT=tm-tt.back();
  MPI::COMM_WORLD.Reduce(&DT,&DT_max,1,MPI_DOUBLE,MPI_MAX,0);
  if(mpirank==0) cout << "tdn elapse time: " << DT_max << " sec" << endl;
  tt.push_back(tm);
  MPI::COMM_WORLD.Barrier();


  //// Fourier transform


  size_t NCs = 2*Ns-1;
  size_t NCsPerCol=NsegPerCol*NCs;

  size_t NCsf = floor((double)NCs/2.0)+1;
  size_t NCsfPerCol=NsegPerCol*NCsf;

  /*
  fftwf_complex ***ssf;
  ssf=new fftwf_complex**[Nr_row];
  for(int ir=0; ir!=Nr_row; ++ir) {
    ssf[ir] = new fftwf_complex*[NsegPerCol];

    for(int iseg=0; iseg!=NsegPerCol; ++iseg) {
      ssf[ir][iseg]= new fftwf_complex[Ns]; //(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*Ns);
    }
  }
  */

  /*
  fftwf_complex **ssf=new fftwf_complex*[Nr_row];
  for(int ir=0; ir!=Nr_row; ++ir)
    ssf[ir]=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*NsegPerCol*Ns);
  */

  fftwf_complex *ssf=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*Nr_row*NCsfPerCol);

  fftwf_plan pfft;
  
  float *in=(float*) fftwf_malloc(sizeof(float)*NCs);
  for(size_t i=0; i!=NCs; ++i) in[i]=0.0;

  pfft=fftwf_plan_dft_r2c_1d(NCs,in,ssf,FFTW_ESTIMATE);

  for(int ir=0; ir!=Nr_row; ++ir) {
    for(int iseg=0; iseg!=NsegPerCol; ++iseg) {
      for(size_t i=0; i!=Ns; ++i) in[i]=ss[ir*NsPerCol+iseg*Ns+i];
      fftwf_execute_dft_r2c(pfft,in,&ssf[ir*NCsfPerCol+iseg*NCsf]);
    }
  }
  fftwf_destroy_plan(pfft);
  fftwf_free(in);
  delete[] ss;

  /*
  if (colrank==12) {
    if (rowrank==0) {
      ofstream fou("psin.ou");
      for(int i=NCs; i!=2*NCs; ++i) {
	fou << creal(ssf[NCsPerCol+i]) << " " << cimag(ssf[NCsPerCol+i]) << endl;
      }
      fou.close();
    }
  }
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Abort(1);
  */


  MPI::COMM_WORLD.Barrier();
  tm=MPI::Wtime();
  DT=tm-tt.back();
  MPI::COMM_WORLD.Reduce(&DT,&DT_max,1,MPI_DOUBLE,MPI_MAX,0);
  if(mpirank==0) cout << "fft elapse time: " << DT_max << " sec" << endl;
  tt.push_back(tm);
  MPI::COMM_WORLD.Barrier();


  //// frequency domain normalization

  const double PI  =3.141592653589793238463;

  int winf = atoi(argv[7]);

  /*
  double f1=(double)atol(argv[8]);
  double f2=(double)atol(argv[9]);
  double f3=(double)atol(argv[10]);
  double f4=(double)atol(argv[11]);

  cout << "f1,2,3,4: " << f1 << " " << f2 << " " << f3 << " " << f4 << " NCsf: " << NCsf << endl;
  */
  
  float* spec = new float[NCsf];
  double* ffc = new double[NCsf];
  double* wfc = new double[NCsf];

  for(int ir=0; ir!=Nr_row; ++ir) {
    for(int iseg=0; iseg!=NsegPerCol; ++iseg) {

      // amplitude spectrum

      for(size_t i=0; i!=NCsf; ++i) 
	spec[i]=sqrt(cabs(ssf[ir*NCsfPerCol+iseg*NCsf+i]*conj(ssf[ir*NCsfPerCol+iseg*NCsf+i])));

      // moving average 

      ffc[0]=spec[0];
      for(size_t i=1; i!=NCsf; ++i)
	ffc[i]=ffc[i-1]+spec[i];

      // weights

      for(size_t i=winf; i!=NCsf-winf; ++i) {
	double tmp=ffc[i+winf]-ffc[i-winf];
	wfc[i]=(tmp==0?1:tmp);	
      }
      for(size_t i=0; i!=winf; ++i)
	wfc[i]=wfc[winf];
      for(size_t i=NCsf-winf; i!=NCsf; ++i)
	wfc[i]=wfc[NCsf-winf-1];

      // normalize

      for(size_t i=0; i!=NCsf; ++i)
	ssf[ir*NCsfPerCol+iseg*NCsf+i] /= wfc[i];

      // filter with taper f1, f2, f3, f4

      /*
      for(size_t i=0; i!=(size_t)f1; ++i)
	ssf[ir*NCsfPerCol+iseg*NCsf+i]=0.0;

      for(size_t i=(size_t)f1; i!=(size_t)f2; ++i) {
	double w=pow(sin((i-f1)/(f2-f1)*PI/2),2);
	ssf[ir*NCsfPerCol+iseg*NCsf+i] *= w;
      }
      
      for(size_t i=(size_t)f3; i!=(size_t)f4; ++i) {
	double w=pow(sin((i-f4)/(f3-f4)*PI/2),2);
	ssf[ir*NCsfPerCol+iseg*NCsf+i] *= w;
      }

      for(size_t i=(size_t)f4; i!=NCsf; ++i)
	ssf[ir*NCsfPerCol+iseg*NCsf+i] = 0.0;
      */
      
    } // for iseg
  } // for ir

  delete[] spec;
  delete[] ffc;
  delete[] wfc;


  //// local xcorr

  MPI::COMM_WORLD.Barrier();
  tm=MPI::Wtime();
  tt.push_back(tm);
  MPI::COMM_WORLD.Barrier();


  int Nx_max=(Nr_row-1)*Nr_row/2;

  fftwf_complex *cssf=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*NCsf);
  float *css=(float*)fftwf_malloc(sizeof(float)*NCs);

  fftwf_plan pifft;
  pifft=fftwf_plan_dft_c2r_1d(NCs,cssf,css,FFTW_ESTIMATE);
  
  //int Nxl = (Nr_row-1)*Nr_row/2;

  /*
  fftwf_complex ***cssf=new fftwf_complex**[Nxl];
  for(int ic=0; ic!=Nxl; ++ic) {
    cssf[ic]=new fftwf_complex*[NsegPerCol];
    for(int iseg=0; iseg!=NsegPerCol; ++iseg)
      cssf[ic][iseg]=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*NCs);
  }
  */

  vector<int> pairlst_col_all;
  vector<int> pairlst;
  float* cssl=new float[Nx_max*NCs]();

  int ic=0;

  for(int ir1=0; ir1!=Nr_row-1; ++ir1) {
    for(int ir2=ir1+1; ir2!=Nr_row; ++ir2) {

      if(mpirank==0)
	if(ic%10 ==0) cout << "ic = " << ic << " out of " << Nx_max << endl;

      pairlst.push_back(rlst_row[ir1]);
      pairlst.push_back(rlst_row[ir2]);

      for(int iseg=0; iseg!=NsegPerCol; ++iseg) {

	for(size_t i=0; i!=NCsf; ++i)
	  cssf[i]=ssf[ir1*NCsfPerCol+iseg*NCsf+i]*conj(ssf[ir2*NCsfPerCol+iseg*NCsf+i]);

	fftwf_execute_dft_c2r(pifft,cssf,css);

	// normalize time-domain css

	float cmax=0.0;
	for(size_t i=0; i!=NCs; ++i) {
	  if(fabs(css[i])>cmax) cmax=fabs(css[i]);
	}
	for(size_t i=0; i!=NCs; ++i) {
	  css[i] /= (cmax==0?1:cmax);
	}
	
	// local stack 

	for(size_t i=0; i!=NCs; ++i)
	  cssl[ic*NCs+i] += css[i];

      } // for iseg

      ++ic;
    } // for ir2
  } // for ir1


  MPI::COMM_WORLD.Barrier();
  tm=MPI::Wtime();
  DT=tm-tt.back();
  MPI::COMM_WORLD.Reduce(&DT,&DT_max,1,MPI_DOUBLE,MPI_MAX,0);
  if(mpirank==0) cout << "local xcorr local stack elapse time: " << DT_max << " sec" << endl;
  tt.push_back(tm);
  MPI::COMM_WORLD.Barrier();


  // global stack

  //int Nxg = ceil((float)Nr*(float)(Nr-1)/2.0/(float)mpisize);
  //float* cssg=new float[Nxg*NCs]();
  vector<float> cssg;

  vector<int> NxPerCol(dims[1],0);
  int icol=0;
  for(size_t i=0; i!=Nx_max; ++i) {
    NxPerCol[icol] += 1;
    /*
    if(mpirank==0) {
      cout << "i: " << i << " icol: " << icol << " NxPerCol: " << NxPerCol[icol] << endl;
    }
    */
    if (icol<dims[1]-1) ++icol;
    else icol=0;
  }
  cssg.resize(NxPerCol[rowrank]*NCs,0.0);

  vector<int> Nxbeg(dims[1],0);
  for(int i=1; i!=dims[1]; ++i)
    Nxbeg[i]=Nxbeg[i-1]+NxPerCol[i-1];
  vector<int>::const_iterator b0=pairlst.begin()+Nxbeg[rowrank]*2;
  vector<int>::const_iterator e0=pairlst.begin()+Nxbeg[rowrank]*2+NxPerCol[rowrank]*2;
  vector<int> pairlst_col(b0,e0);
  pairlst_col_all.insert(pairlst_col_all.end(),pairlst_col.begin(),pairlst_col.end());

  /*
  if(rowrank==0) {
    cout << "colrank: " << colrank << " ";
    for(int i=0; i!=dims[1]; ++i) cout << NxPerCol[i] << " " ;
    cout << "pairlst_col_all.size: " << pairlst_col_all.size() << endl;
  }
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Abort(1);
  */

  vector<size_t> NxB(dims[1],0);
  for (int i=1; i!=dims[1]; ++i)
    NxB[i]=NxB[i-1]+NxPerCol[i-1]*NCs;
  
  MPI::COMM_WORLD.Barrier();
  double tstack0=MPI::Wtime();

  for(int i=0; i!=dims[1]; ++i) {
    comrow.Reduce(&cssl[NxB[i]],&cssg[0],NxPerCol[i]*NCs,MPI_FLOAT,MPI_SUM,i);
  }

  double tstack1=MPI::Wtime();
  double dtstack=tstack1-tstack0;
  double dtstack_max=0.0;
  MPI::COMM_WORLD.Reduce(&dtstack,&dtstack_max,1,MPI_DOUBLE,MPI_MAX,0);
  if(mpirank==0) cout << "local xcorr global reduce elapse time: " << dtstack_max << " sec" << endl;
  MPI::COMM_WORLD.Barrier();

  vector<size_t> NxgE(dims[1],0);
  for(int i=0; i!=dims[1]; ++i) 
    NxgE[i]=NxPerCol[i]*NCs;

  delete[] cssl;

  /*
  if(rowrank==0 && colrank==0){
    ofstream fou("psin.ou00", ios::out|ios::binary);
    fou.write(reinterpret_cast<char*>(cssg),sizeof(float)*(NCs));
    fou.close();
  }
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Abort(1);
  */

  MPI::COMM_WORLD.Barrier();
  tm=MPI::Wtime();
  DT=tm-tt.back();
  MPI::COMM_WORLD.Reduce(&DT,&DT_max,1,MPI_DOUBLE,MPI_MAX,0);
  if(mpirank==0) cout << "local xcorr global stack elapse time: " << DT_max << " sec" << endl;
  tt.push_back(tm);
  MPI::COMM_WORLD.Barrier();


  int Nx_all=Nr*(Nr-1)/2;
  int Nx_done=0;
  int Nx_done_all=0;
  Nx_done=pairlst_col_all.size()/2;  
  com2d.Allreduce(&Nx_done,&Nx_done_all,1,MPI_INT,MPI_SUM);
  if(mpirank==0) {
    cout << "completed "<< Nx_done_all<<" out of " << Nx_all << " xcorr" << endl;
    cout << "..................................................." << endl;
  }
  
  /*
  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Abort(1);
  */


  //// neighbor xcorr



  MPI::COMM_WORLD.Barrier();
  tm=MPI::Wtime();
  tt.push_back(tm);
  MPI::COMM_WORLD.Barrier();

  int Nr_row_max, Nr_row1=Nr_row;
  fftwf_complex *ssf1;
  vector<int> rlst_row1;
  int up, down;
  int Nshift=0;
  bool ssf1_alloc=false;
  
  //cout << "check 1 rank: " << mpirank << " " << Nx_done_all << " " << Nx_all << endl;
  if (Nx_done_all < Nx_all) {

    ssf1_alloc=true;
    
    MPI::COMM_WORLD.Allreduce(&Nr_row,&Nr_row_max,1,MPI_INT,MPI_MAX);
    //cout << "check 1.1 rank: " << mpirank << " " << Nr_row_max << endl;
    rlst_row1.resize(Nr_row_max,0);
    
    ssf1=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*Nr_row_max*NCsfPerCol);  
    for(size_t i=0; i!=Nr_row*NCsfPerCol; ++i) ssf1[i]=ssf[i];
    
    for(int i=0; i!=Nr_row; ++i) rlst_row1[i]=rlst_row[i];
    
    com2d.Shift(0,1,down,up);

  }


  //for(int ishift=0; ishift!=Nshift; ++ishift) {
  //do {
  //cout << "check 2 rank: " << mpirank << " " << Nx_done_all << " " << Nx_all << endl;
  while (Nx_done_all<Nx_all) {
    
    MPI::COMM_WORLD.Barrier();
    tm=MPI::Wtime();
    tt.push_back(tm);

    com2d.Sendrecv_replace(&Nr_row1,1,MPI_INT,up,0,down,0);
    com2d.Sendrecv_replace(ssf1,Nr_row_max*NCsfPerCol,MPI_C_FLOAT_COMPLEX,up,0,down,0);
    com2d.Sendrecv_replace(&rlst_row1[0],Nr_row_max,MPI_INT,up,0,down,0);

    MPI::COMM_WORLD.Barrier();
    tm=MPI::Wtime();
    DT=tm-tt.back();
    MPI::COMM_WORLD.Reduce(&DT,&DT_max,1,MPI_DOUBLE,MPI_MAX,0);
    if(mpirank==0) cout << "SendRecv elapse time: " << DT_max << " sec" << endl;
    tt.push_back(tm);
    MPI::COMM_WORLD.Barrier();

    /*
    if(mpirank==0) cout << "shift 1: " << Nr_row << " " << Nr_row1 << endl;
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Abort(1);
    */
    /*
    if (colrank==13) {
      if (rowrank==0) {
	ofstream fou("psin_fft1.ou");
	for(int i=NCs; i!=2*NCs; ++i) {
	  fou << creal(ssf1[NCsPerCol+i]) << " " << cimag(ssf1[NCsPerCol+i]) << endl;
	}
	fou.close();
      }
    }
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Abort(1);
    */

    //double tb=MPI::Wtime();

    Nx_max=Nr_row*Nr_row1;
    cssl=new float[Nx_max*NCs]();
    pairlst.clear();

    ic = 0;
    for(int ir1=0; ir1!=Nr_row; ++ir1) {
      for(int ir2=0; ir2!=Nr_row1; ++ir2) {

	if(mpirank==0) {
	  if(ic%10 ==0) cout << "ic = " << ic << " out of " << Nx_max << endl;
	}

	pairlst.push_back(rlst_row[ir1]);
	pairlst.push_back(rlst_row1[ir2]);

	for(int iseg=0; iseg!=NsegPerCol; ++iseg) {
	  
	  for(size_t i=0; i!=NCsf; ++i) {
	    cssf[i]=ssf[ir1*NCsfPerCol+iseg*NCsf+i]*conj(ssf1[ir2*NCsfPerCol+iseg*NCsf+i]);
	  }
	  fftwf_execute_dft_c2r(pifft,cssf,css);
	  
	  // normalize time-domain css
	  
	  float cmax=0.0;
	  for(size_t i=0; i!=NCs; ++i) {
	    if(fabs(css[i])>cmax) cmax=fabs(css[i]);
	  }
	  for(size_t i=0; i!=NCs; ++i) {
	    css[i] /= (cmax==0?1:cmax);
	  }
	  
	  // local stack 
	  
	  for(size_t i=0; i!=NCs; ++i)
	    cssl[ic*NCs+i] += css[i];
	  
	} // for iseg

	++ic;
      } // for ir2
    } // for ir1

    /*
    double te=MPI::Wtime();
    double dteb=te-tb;
    cout << "rank: " << mpirank << " dteb= " << dteb << endl;
    */

    MPI::COMM_WORLD.Barrier();
    tm=MPI::Wtime();
    DT=tm-tt.back();
    MPI::COMM_WORLD.Reduce(&DT,&DT_max,1,MPI_DOUBLE,MPI_MAX,0);
    if(mpirank==0) cout << "neighbor xcorr local stack elapse time: " << DT_max << " sec" << endl;
    tt.push_back(tm);
    MPI::COMM_WORLD.Barrier();

    
    // global stack
    

    //double tb=MPI::Wtime();

    NxPerCol.assign(dims[1],0);
    icol=0;
    for(size_t i=0; i!=Nx_max; ++i) {
      NxPerCol[icol] += 1;
      /*
      if(mpirank==0) {
	cout << "i: " << i << " icol: " << icol << " NxPerCol: " << NxPerCol[icol] << endl;
      }
      */
      if (icol<dims[1]-1) ++icol;
      else icol=0;
    }
    /*
    if(mpirank==0) {
      cout << "Nx_max: " << Nx_max << endl;
      for(int i=0; i!=dims[1]; ++i)
	cout << NxPerCol[i] << " " ;
      cout << endl;
    }
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Abort(1);
    */
    cssg.resize(cssg.size()+NxPerCol[rowrank]*NCs);
    
    Nxbeg.assign(dims[1],0);
    for(int i=1; i!=dims[1]; ++i)
      Nxbeg[i]=Nxbeg[i-1]+NxPerCol[i-1];
    b0=pairlst.begin()+Nxbeg[rowrank]*2;
    e0=pairlst.begin()+Nxbeg[rowrank]*2+NxPerCol[rowrank]*2;
    pairlst_col.assign(b0,e0);
    pairlst_col_all.insert(pairlst_col_all.end(),pairlst_col.begin(),pairlst_col.end());

    /*
    if(colrank==0) {
      cout << "rowrank: " << rowrank << " NxPerCol: " << NxPerCol[rowrank] <<  " size: " << pairlst_col_all.size() << " ";
      for(int i=0; i!=pairlst_col_all.size(); ++i) cout << pairlst_col_all[i] << " ";
      cout << endl;
    }
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Abort(1);
    */

    /*
    if(mpirank==0) {
      for(int i=0; i!=dims[1]; ++i) cout << NxPerCol[i] << " " ;
      cout << endl;
    }
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Abort(1);
    */

    
    NxB.assign(dims[1],0);
    for (int i=1; i!=dims[1]; ++i)
      NxB[i]=NxB[i-1]+NxPerCol[i-1]*NCs;

    /*
    double te=MPI::Wtime();
    double dteb=te-tb;
    cout << "rank: " << mpirank << " dteb= " << dteb << endl;
    */
    

    MPI::COMM_WORLD.Barrier();
    double tstack2=MPI::Wtime();

    for(int i=0; i!=dims[1]; ++i) {
      comrow.Reduce(&cssl[NxB[i]],&cssg[NxgE[i]],NxPerCol[i]*NCs,MPI_FLOAT,MPI_SUM,i);
    }

    double tstack3=MPI::Wtime();
    double dtstack1=tstack3-tstack2;
    double dtstack1_max=0;
    MPI::COMM_WORLD.Reduce(&dtstack1,&dtstack1_max,1,MPI_DOUBLE,MPI_MAX,0);
    if(mpirank==0) cout << "neighbor xcorr global reduce elapse time: " << dtstack1_max << " sec" << endl;
    MPI::COMM_WORLD.Barrier();

    tm=MPI::Wtime();
    DT=tm-tt.back();
    MPI::COMM_WORLD.Reduce(&DT,&DT_max,1,MPI_DOUBLE,MPI_MAX,0);
    if(mpirank==0) cout << "neighbor xcorr global stack elapse time: " << DT_max << " sec" << endl;
    tt.push_back(tm);
    MPI::COMM_WORLD.Barrier();

    for(int i=0; i!=dims[1]; ++i) 
      NxgE[i]+=NxPerCol[i]*NCs;
    
    delete[] cssl;

    ++Nshift;    

    /*
    cout << "rowrank: " << rowrank << " colrank: " << colrank << " 2 * xcorr #: " << pairlst_col_all.size() << endl;
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Abort(1);
    */

    Nx_done=pairlst_col_all.size()/2;
    com2d.Allreduce(&Nx_done,&Nx_done_all,1,MPI_INT,MPI_SUM);

    if(mpirank==0) {
      cout << "completed "<< Nx_done_all<<" out of " << Nx_all << " xcorr" << endl;
      cout << "shuffle number: " << Nshift << " completed..." << endl;
      cout << "..................................................." << endl;
    }

    MPI::COMM_WORLD.Barrier();

    /*
    if (Nshift==2) { // for weak-scaling figure
      MPI::COMM_WORLD.Barrier();
      MPI::COMM_WORLD.Abort(1);
    }
    */

  }  // while for ishift
  
  fftwf_destroy_plan(pifft);
  fftwf_free(css);
  fftwf_free(cssf);
  fftwf_free(ssf);
  if(ssf1_alloc) fftwf_free(ssf1);
  //cout << "check 3 rank: " << mpirank << " " << Nx_done_all << " " << Nx_all << endl;
  
  // for scaling tests, comment out the next 2 lines for real run
  //MPI::COMM_WORLD.Barrier();
  //MPI::COMM_WORLD.Abort(1);


  //// write out cssg

  tm=MPI::Wtime();
  tt.push_back(tm);


  size_t szout=cssg.size();
  if(szout>0) {

    string outpth(argv[8]); //12]);

    stringstream iss;
    iss << setw(2) << setfill('0') << rowrank;
    string colid; iss >> colid;
    iss.str(string()); iss.clear();
    iss << setw(2) << setfill('0') << colrank;
    string rowid; iss >> rowid;

    ofstream fou((outpth+"/psin.ou"+rowid+colid).c_str(), ios::out|ios::binary);
    fou.write(reinterpret_cast<char*>(&cssg[0]), szout*sizeof(float));
    fou.close();

    int psz=pairlst_col_all.size();
    ofstream fou1((outpth+"/psin.idx"+rowid+colid).c_str(),ios::out|ios::binary);
    fou1.write(reinterpret_cast<char*>(&pairlst_col_all[0]),psz*sizeof(int));
    fou1.close();

  }

  //delete[] cssg;

  tm=MPI::Wtime();
  DT=tm-tt.back();
  MPI::COMM_WORLD.Reduce(&DT,&DT_max,1,MPI_DOUBLE,MPI_MAX,0);
  if(mpirank==0) cout << "write output elapse time: " << DT_max << " sec" << endl;
  tt.push_back(tm);

  // write out timing array

  if(mpirank==0) {
    ofstream tfou("timing.dat");
    for(int i=0; i!=tt.size(); ++i)
      tfou << scientific << setprecision(20) << tt[i] << endl;
    tfou.close();
  }


  MPI::Finalize();
  return 0;
}

    
