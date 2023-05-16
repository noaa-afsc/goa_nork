#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
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
#include <ageage.htp>

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
  nobs.allocate("nobs");
  age.allocate(1,nobs,"age");
  ape.allocate(1,nobs,"ape");
  n.allocate(1,nobs,"n");
}

void model_parameters::initializationfunction(void)
{
  sigma1.set_initial_value(0.5);
  sigma2.set_initial_value(5);
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
  sigma1.allocate(1,"sigma1");
  sigma2.allocate(2,"sigma2");
  Perc_Corr.allocate(1,nobs,"Perc_Corr");
  #ifndef NO_AD_INITIALIZE
    Perc_Corr.initialize();
  #endif
  Perc_Corr1.allocate(1,nobs,"Perc_Corr1");
  #ifndef NO_AD_INITIALIZE
    Perc_Corr1.initialize();
  #endif
  Perc_Corr2.allocate(1,nobs,"Perc_Corr2");
  #ifndef NO_AD_INITIALIZE
    Perc_Corr2.initialize();
  #endif
  Phat.allocate(1,nobs,"Phat");
  #ifndef NO_AD_INITIALIZE
    Phat.initialize();
  #endif
  sigma_a.allocate(1,nobs,"sigma_a");
  sigma_inc.allocate("sigma_inc");
  #ifndef NO_AD_INITIALIZE
  sigma_inc.initialize();
  #endif
  RSS.allocate(1,nobs,"RSS");
  #ifndef NO_AD_INITIALIZE
    RSS.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::userfunction(void)
{
  f =0.0;
  get_A_SD_est();
  evaluate_the_objective_function();
}

void model_parameters::get_A_SD_est(void)
{
  sigma_inc = (sigma2-sigma1)/(age(nobs)-age(1));
  sigma_a(1) = sigma1;
  for (int i=2;i<=nobs;i++){
   sigma_a(i) = sigma1+(age(i)-age(1))*sigma_inc;}
  for (int i=1;i<=nobs;i++){
  Perc_Corr(i) = cumd_norm(0.5/sigma_a(i))-cumd_norm((-0.5)/sigma_a(i));
  Perc_Corr1(i) = cumd_norm((-0.5)/sigma_a(i))-cumd_norm((-1.5)/sigma_a(i));
  Perc_Corr2(i) = cumd_norm((-1.5)/sigma_a(i))-cumd_norm((-2.5)/sigma_a(i));}
  Phat = square(Perc_Corr)+2*square(Perc_Corr1)+2*square(Perc_Corr2);
}

void model_parameters::set_runtime(void)
{
  dvector temp("{1e-7,1e-7}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::evaluate_the_objective_function(void)
{
  RSS.initialize();
  for (int i=1;i<=nobs;i++){
   if(n(i)>0){
   RSS(i) = pow(n(i),0.5)*square(Phat(i)-ape(i));}}
  f = sum(RSS);
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

void model_parameters::report(const dvector& gradients){}

void model_parameters::final_calcs(void){}

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
    gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
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
