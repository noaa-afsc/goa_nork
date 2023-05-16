#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  #include <admodel.h>
  #undef REPORT
  #define write_R(object) mysum << #object "\n" << object << endl;
  ofstream mysum("Est_out.rep");
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <vbl.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
  nages.allocate("nages");
  ages.allocate(1,nages,"ages");
  Lbar_obs.allocate(1,nages,"Lbar_obs");
  SD_Lbar.allocate(1,nages,"SD_Lbar");
}

void model_parameters::initializationfunction(void)
{
  Linf.set_initial_value(50);
  k.set_initial_value(0.1);
  t0.set_initial_value(0.1);
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  Linf.allocate(1,"Linf");
  k.allocate(2,"k");
  t0.allocate(3,"t0");
  Lbar_est.allocate(1,nages,"Lbar_est");
  #ifndef NO_AD_INITIALIZE
    Lbar_est.initialize();
  #endif
  yvar.allocate(1,nages,"yvar");
  #ifndef NO_AD_INITIALIZE
    yvar.initialize();
  #endif
  yconst.allocate(1,nages,"yconst");
  #ifndef NO_AD_INITIALIZE
    yconst.initialize();
  #endif
  RSS.allocate(1,nages,"RSS");
  #ifndef NO_AD_INITIALIZE
    RSS.initialize();
  #endif
  jnll.allocate("jnll");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::userfunction(void)
{
  jnll =0.0;
  get_Lbar_est();
  evaluate_the_objective_function();
}

void model_parameters::get_Lbar_est(void)
{
  for (int i=1;i<=nages;i++){
   Lbar_est(i) = Linf*(1-exp(-1.0*k*(ages(i)-t0)));
   yvar(i) = log(1.+square(SD_Lbar(i))/square(Lbar_obs(i)));
   yconst(i) = log(2.0*M_PI*yvar(i)*square(Lbar_obs(i)));}
}

void model_parameters::evaluate_the_objective_function(void)
{
  for (int i=1;i<=nages;i++){
   RSS(i) = 0.5*(yconst(i)+square(log(Lbar_est(i))-log(Lbar_obs(i)))/yvar(i));}
  jnll = sum(RSS);
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report << "Linf "<<Linf<<endl;
  report << "k "<<k<<endl;
  report << "t0 "<<t0<<endl;
  report << "Obj_fun "<<jnll<<endl;
}

void model_parameters::final_calcs()
{
  write_R(ages);
  write_R(Lbar_obs);
  write_R(Lbar_est);
  mysum.close();
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
