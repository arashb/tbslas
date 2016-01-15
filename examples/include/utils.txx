/**
 * \author Dhairya Malhotra, dhairya.malhotra@gmail.com
 * \author Arash Bakhtiari
 */

#include <utils/metadata.h>
#include <utils/sim_config.h>

typedef tbslas::MetaData<std::string,
                         std::string,
                         std::string> MetaData_t;

template <class FMM_Mat_t>
void CheckFMMOutput(pvfmm::FMM_Tree<FMM_Mat_t>* mytree, const pvfmm::Kernel<typename FMM_Mat_t::Real_t>* mykernel){
  if(mykernel==NULL) return;

  // Find out number of OMP thereads.
  int np=omp_get_max_threads();

  // Find out my identity in the default communicator
  int myrank, p;
  MPI_Comm c1=MPI_COMM_WORLD;
  MPI_Comm_rank(c1, &myrank);
  MPI_Comm_size(c1,&p);

  typedef typename FMM_Mat_t::FMMData FMM_Data_t;
  typedef typename FMM_Mat_t::FMMNode_t FMMNode_t;
  typedef typename FMM_Mat_t::Real_t Real_t;

  // Read source data.
  std::vector<Real_t> src_coord;
  std::vector<Real_t> src_value;
  FMMNode_t* n=static_cast<FMMNode_t*>(mytree->PreorderFirst());
  while(n!=NULL){
    if(n->IsLeaf() && !n->IsGhost()){
      pvfmm::Vector<Real_t>& coord_vec=n->src_coord;
      pvfmm::Vector<Real_t>& value_vec=n->src_value;
      for(size_t i=0;i<coord_vec.Dim();i++) src_coord.push_back(coord_vec[i]);
      for(size_t i=0;i<value_vec.Dim();i++) src_value.push_back(value_vec[i]);
    }
    n=static_cast<FMMNode_t*>(mytree->PreorderNxt(n));
  }
  long long glb_src_cnt=0, src_cnt=src_coord.size()/3;
  MPI_Allreduce(&src_cnt, &glb_src_cnt, 1, MPI_LONG_LONG, MPI_SUM, c1);
  long long glb_val_cnt=0, val_cnt=src_value.size();
  MPI_Allreduce(&val_cnt, &glb_val_cnt, 1, MPI_LONG_LONG, MPI_SUM, c1);
  if(glb_src_cnt==0) return;

  int dof=glb_val_cnt/glb_src_cnt/mykernel->ker_dim[0];
  int trg_dof=dof*mykernel->ker_dim[1];

  // Read target data.
  std::vector<Real_t> trg_coord;
  std::vector<Real_t> trg_poten_fmm;
  long long trg_iter=0;
  size_t step_size=1+glb_src_cnt*glb_src_cnt*1e-9/p;
  n=static_cast<FMMNode_t*>(mytree->PreorderFirst());
  while(n!=NULL){
    if(n->IsLeaf() && !n->IsGhost()){
      pvfmm::Vector<Real_t>& coord_vec=n->trg_coord;
      pvfmm::Vector<Real_t>& poten_vec=n->trg_value;
      for(size_t i=0;i<coord_vec.Dim()/3          ;i++){
        if(trg_iter%step_size==0){
          for(int j=0;j<3        ;j++) trg_coord    .push_back(coord_vec[i*3        +j]);
          for(int j=0;j<trg_dof  ;j++) trg_poten_fmm.push_back(poten_vec[i*trg_dof  +j]);
        }
        trg_iter++;
      }
    }
    n=static_cast<FMMNode_t*>(mytree->PreorderNxt(n));
  }
  int trg_cnt=trg_coord.size()/3;
  int send_cnt=trg_cnt*3;
  std::vector<int> recv_cnts(p), recv_disp(p,0);
  MPI_Allgather(&send_cnt    , 1, MPI_INT,
                &recv_cnts[0], 1, MPI_INT, c1);
  pvfmm::omp_par::scan(&recv_cnts[0], &recv_disp[0], p);
  int glb_trg_cnt=(recv_disp[p-1]+recv_cnts[p-1])/3;
  std::vector<Real_t> glb_trg_coord(glb_trg_cnt*3);
  MPI_Allgatherv(&trg_coord[0]    , send_cnt                    , pvfmm::par::Mpi_datatype<Real_t>::value(),
                 &glb_trg_coord[0], &recv_cnts[0], &recv_disp[0], pvfmm::par::Mpi_datatype<Real_t>::value(), c1);
  if(glb_trg_cnt==0) return;

  //Direct N-Body.
  std::vector<Real_t> trg_poten_dir(glb_trg_cnt*trg_dof ,0);
  std::vector<Real_t> glb_trg_poten_dir(glb_trg_cnt*trg_dof ,0);
  pvfmm::Profile::Tic("N-Body Direct",&c1,false,1);
  #pragma omp parallel for
  for(int i=0;i<np;i++){
    size_t a=(i*glb_trg_cnt)/np;
    size_t b=((i+1)*glb_trg_cnt)/np;
    mykernel->ker_poten(&src_coord[0], src_cnt, &src_value[0], dof, &glb_trg_coord[a*3], b-a, &trg_poten_dir[a*trg_dof  ],NULL);
  }
  MPI_Allreduce(&trg_poten_dir[0], &glb_trg_poten_dir[0], trg_poten_dir.size(), pvfmm::par::Mpi_datatype<Real_t>::value(), pvfmm::par::Mpi_datatype<Real_t>::sum(), c1);
  pvfmm::Profile::Toc();

  //Compute error.
  {
    Real_t max_=0;
    Real_t max_err=0;
    for(size_t i=0;i<trg_poten_fmm.size();i++){
      Real_t err=fabs(glb_trg_poten_dir[i+(recv_disp[myrank]/3)*trg_dof]-trg_poten_fmm[i]);
      Real_t max=fabs(glb_trg_poten_dir[i+(recv_disp[myrank]/3)*trg_dof]);
      if(err>max_err) max_err=err;
      if(max>max_) max_=max;
    }
    Real_t glb_max, glb_max_err;
    MPI_Reduce(&max_   , &glb_max    , 1, pvfmm::par::Mpi_datatype<Real_t>::value(), pvfmm::par::Mpi_datatype<Real_t>::max(), 0, c1);
    MPI_Reduce(&max_err, &glb_max_err, 1, pvfmm::par::Mpi_datatype<Real_t>::value(), pvfmm::par::Mpi_datatype<Real_t>::max(), 0, c1);
    if(!myrank){
      std::cout<<"Maximum Absolute Error ["<<mykernel->ker_name<<"] :  "<<std::scientific<<glb_max_err<<'\n';
      std::cout<<"Maximum Relative Error ["<<mykernel->ker_name<<"] :  "<<std::scientific<<glb_max_err/glb_max<<'\n';
    }
  }
}


template <class FMMTree_t,
          class Functor_t>
void CheckChebOutput(FMMTree_t* mytree,
                     Functor_t fn_poten,
                     int fn_dof,
                     typename FMMTree_t::Real_t& al2,
                     typename FMMTree_t::Real_t& rl2,
                     typename FMMTree_t::Real_t& ali,
                     typename FMMTree_t::Real_t& rli,
                     std::string t_name) {
  typedef typename FMMTree_t::Node_t FMMNode_t;
  typedef typename FMMTree_t::Real_t Real_t;

  MPI_Comm c1=*mytree->Comm();
  pvfmm::Profile::Tic((std::string("Compute Error ")+t_name).c_str(),&c1,true,1);

  int myrank; MPI_Comm_rank(c1, &myrank);
  FMMNode_t* r_node=static_cast<FMMNode_t*>(mytree->RootNode());
  int dof=r_node->DataDOF()/fn_dof;
  std::vector<FMMNode_t*> nodes;
  {
    FMMNode_t* node=static_cast<FMMNode_t*>(mytree->PreorderFirst());
    while(node!=NULL){
      if(node->IsLeaf() && !node->IsGhost()) nodes.push_back(node);
      node=static_cast<FMMNode_t*>(mytree->PreorderNxt(node));
    }
    if(nodes.size()==0) return;
  }

  int lcl_nodes_size = nodes.size();
  int glb_nodes_size = 0;
  MPI_Allreduce(&lcl_nodes_size, &glb_nodes_size, 1,
                MPI_INT, MPI_SUM, c1);
  int cheb_deg=nodes[0]->ChebDeg();
  std::vector<Real_t> cheb_nds=pvfmm::cheb_nodes<Real_t>(cheb_deg/3+1, 1);
  for(size_t i=0;i<cheb_nds.size();i++) cheb_nds[i]=2.0*cheb_nds[i]-1.0;
  std::vector<Real_t> cheb_pts=pvfmm::cheb_nodes<Real_t>(cheb_deg/3+1, COORD_DIM);
  int n_pts=cheb_pts.size()/COORD_DIM;
  int omp_p=omp_get_max_threads();

  std::vector<Real_t> glb_err_avg(dof*fn_dof,0);
  { // Determine glb_err_avg.
    std::vector<Real_t> err_avg(omp_p*dof*fn_dof,0);
    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      pvfmm::Vector<Real_t> out; out.SetZero();
      pvfmm::Vector<Real_t> fn_out(dof*fn_dof);
      for(size_t i=(nodes.size()*tid)/omp_p;i<(nodes.size()*(tid+1))/omp_p;i++){
        pvfmm::Vector<Real_t>& cheb_coeff=nodes[i]->ChebData();
        cheb_eval(cheb_coeff, cheb_deg, cheb_nds, cheb_nds, cheb_nds, out);

        Real_t* c=nodes[i]->Coord();
        Real_t s=pow(2.0,-nodes[i]->Depth());
        Real_t s3=s*s*s;
        for(size_t j=0;j<n_pts;j++){
          Real_t coord[3]={c[0]+s*cheb_pts[j*COORD_DIM+0],
                           c[1]+s*cheb_pts[j*COORD_DIM+1],
                           c[2]+s*cheb_pts[j*COORD_DIM+2]};
          fn_out.SetZero();
          for(size_t k=0;k<dof;k++)
            fn_poten(coord,1,&fn_out[k*fn_dof]);
          for(size_t k=0;k<dof*fn_dof;k++){
            Real_t err=out[n_pts*k+j]-fn_out[k];
            err_avg[tid*dof*fn_dof+k]+=err*s3;
          }
        }
      }
      for(size_t k=0;k<dof*fn_dof;k++)
        err_avg[tid*dof*fn_dof+k]/=n_pts;
    }
    for(size_t tid=1;tid<omp_p;tid++)
      for(size_t k=0;k<dof*fn_dof;k++)
        err_avg[k]+=err_avg[tid*dof*fn_dof+k];
    MPI_Allreduce(&err_avg[0], &glb_err_avg[0], dof*fn_dof, pvfmm::par::Mpi_datatype<Real_t>::value(), pvfmm::par::Mpi_datatype<Real_t>::sum(), c1);
  }
  { // Write error to file.
    int nn_x=1;
    int nn_y=201;
    int nn_z=201;
    std::vector<Real_t> x(nn_x,1.0/2.0);
    std::vector<Real_t> y(nn_y);
    std::vector<Real_t> z(nn_z);
    //for(int i=0;i<nn_x;i++)
    //  x[i]=((Real_t)i)/(nn_x-1);
    for(int i=0;i<nn_y;i++)
      y[i]=((Real_t)i)/(nn_y-1);
    for(int i=0;i<nn_z;i++)
      z[i]=((Real_t)i)/(nn_z-1);

    FMMNode_t* r_node=static_cast<FMMNode_t*>(mytree->RootNode());
    int fn_dof=r_node->DataDOF()/dof;

    Real_t* fn_out=new Real_t[dof*fn_dof];
    pvfmm::Matrix<Real_t> M_out    (nn_z*fn_dof*dof,nn_y*nn_x,NULL,true); M_out    .SetZero();
    pvfmm::Matrix<Real_t> M_out_err(nn_z*fn_dof*dof,nn_y*nn_x,NULL,true); M_out_err.SetZero();
    r_node->ReadVal(x,y,z,M_out[0],false);
    for(int l0=0;l0<dof;l0++)
    for(int l1=0;l1<fn_dof;l1++)
    for(int i=0;i<nn_x;i++)
    for(int j=0;j<nn_y;j++)
    for(int k=0;k<nn_z;k++){
      int l=l0*fn_dof+l1;
      Real_t ch_coord[3]={x[i],y[j],z[k]};
      Real_t cheb_val=M_out[k+l*nn_z][i+j*nn_x];
      if(fabs(cheb_val)>1.0e-20){
        fn_poten(ch_coord,1,fn_out);
        Real_t err=(cheb_val-fn_out[l])-glb_err_avg[l];
        M_out_err[k+l*nn_z][i+j*nn_x]=err;
      }
    }
    delete[] fn_out;
    pvfmm::Matrix<Real_t> M_global    (nn_z*fn_dof*dof,nn_y*nn_x,NULL,true);
    pvfmm::Matrix<Real_t> M_global_err(nn_z*fn_dof*dof,nn_y*nn_x,NULL,true);
    MPI_Reduce(&M_out    [0][0], &M_global    [0][0], nn_x*nn_y*nn_z*fn_dof*dof, pvfmm::par::Mpi_datatype<Real_t>::value(), pvfmm::par::Mpi_datatype<Real_t>::sum(), 0, c1);
    MPI_Reduce(&M_out_err[0][0], &M_global_err[0][0], nn_x*nn_y*nn_z*fn_dof*dof, pvfmm::par::Mpi_datatype<Real_t>::value(), pvfmm::par::Mpi_datatype<Real_t>::sum(), 0, c1);

    //std::string fname;
    //fname=std::string("result/"    )+t_name+std::string(".mat");
    //if(!myrank) M_global    .Write(fname.c_str());
    //fname=std::string("result/err_")+t_name+std::string(".mat");
    //if(!myrank) M_global_err.Write(fname.c_str());
  }
  std::vector<Real_t> max    (omp_p,0), l2    (omp_p,0);
  std::vector<Real_t> max_err(omp_p,0), l2_err(omp_p,0);
  #pragma omp parallel for
  for(size_t tid=0;tid<omp_p;tid++){
    pvfmm::Vector<Real_t> out; out.SetZero();
    pvfmm::Vector<Real_t> fn_out(dof*fn_dof);
    for(size_t i=(nodes.size()*tid)/omp_p;i<(nodes.size()*(tid+1))/omp_p;i++){
      pvfmm::Vector<Real_t>& cheb_coeff=nodes[i]->ChebData();
      cheb_eval(cheb_coeff, cheb_deg, cheb_nds, cheb_nds, cheb_nds, out);

      Real_t* c=nodes[i]->Coord();
      Real_t s=pow(2.0,-nodes[i]->Depth());
      Real_t s3=s*s*s;
      for(size_t j=0;j<n_pts;j++){
        Real_t coord[3]={c[0]+s*cheb_pts[j*COORD_DIM+0],
                         c[1]+s*cheb_pts[j*COORD_DIM+1],
                         c[2]+s*cheb_pts[j*COORD_DIM+2]};
        fn_out.SetZero();
        for(size_t k=0;k<dof;k++)
          fn_poten(coord,1,&fn_out[k*fn_dof]);
        for(size_t k=0;k<dof*fn_dof;k++){
          Real_t err=out[n_pts*k+j]-fn_out[k]-glb_err_avg[k];
          if(fabs(fn_out[k])>max    [tid]) max    [tid]=fabs(fn_out[k]);
          if(fabs(err      )>max_err[tid]) max_err[tid]=fabs(err      );
          l2[tid]+=fn_out[k]*fn_out[k]*s3;
          l2_err[tid]+=err*err*s3;
        }
      }
    }
    l2    [tid]/=n_pts;
    l2_err[tid]/=n_pts;
  }
  for(size_t tid=1;tid<omp_p;tid++){
    if(max    [tid]>max    [0]) max    [0]=max    [tid];
    if(max_err[tid]>max_err[0]) max_err[0]=max_err[tid];
    l2    [0]+=l2    [tid];
    l2_err[0]+=l2_err[tid];
  }

  Real_t global_l2, global_l2_err;
  Real_t global_max, global_max_err;
  MPI_Reduce(&l2     [0], &global_l2     , 1, pvfmm::par::Mpi_datatype<Real_t>::value(), pvfmm::par::Mpi_datatype<Real_t>::sum(), 0, c1);
  MPI_Reduce(&l2_err [0], &global_l2_err , 1, pvfmm::par::Mpi_datatype<Real_t>::value(), pvfmm::par::Mpi_datatype<Real_t>::sum(), 0, c1);
  MPI_Reduce(&max    [0], &global_max    , 1, pvfmm::par::Mpi_datatype<Real_t>::value(), pvfmm::par::Mpi_datatype<Real_t>::max(), 0, c1);
  MPI_Reduce(&max_err[0], &global_max_err, 1, pvfmm::par::Mpi_datatype<Real_t>::value(), pvfmm::par::Mpi_datatype<Real_t>::max(), 0, c1);
  al2 = sqrt(global_l2_err);
  rl2 = sqrt(global_l2_err/global_l2);
  ali = global_max_err;
  rli = global_max_err/global_max;

  if(!myrank){
    std::cout<<"Absolute L2 Error ["<<t_name<<"]     :  "<<std::scientific<<sqrt(global_l2_err)<<'\n';
    std::cout<<"Relative L2 Error ["<<t_name<<"]     :  "<<std::scientific<<sqrt(global_l2_err/global_l2)<<'\n';
    std::cout<<"Maximum Absolute Error ["<<t_name<<"]:  "<<std::scientific<<global_max_err<<'\n';
    std::cout<<"Maximum Relative Error ["<<t_name<<"]:  "<<std::scientific<<global_max_err/global_max<<'\n';
  }
  pvfmm::Profile::Toc();
}

void commandline_option_start(int argc, char** argv, const char* help_text){
  char help[]="--help";
  for(int i=0;i<argc;i++){
    if(!strcmp(argv[i],help)){
      if(help_text!=NULL) std::cout<<help_text<<'\n';
      std::cout<<"Usage:\n";
      std::cout<<"  "<<argv[0]<<" [options]\n\n";
    }
  }
}

const char* commandline_option(int argc, char** argv, const char* opt, const char* def_val, bool required, const char* err_msg){
  char help[]="--help";
  for(int i=0;i<argc;i++){
    if(!strcmp(argv[i],help)){
      std::cout<<"        "<<err_msg<<'\n';
      return def_val;
    }
  }
  for(int i=0; i<argc; i++) {
    if(!strcmp(argv[i], opt)) {
      MetaData_t::AddMetaData(opt, argv[(i+1)%argc], err_msg);
      return argv[(i+1)%argc];
    }
  }
  if(required){
    std::cout<<"Missing: required option\n"<<"    "<<err_msg<<"\n\n";
    std::cout<<"To see usage options\n"<<"    "<<argv[0]<<" --help\n\n";
    exit(0);
  }

  std::string val("0");
  if(def_val)
    val = std::string(def_val);

  MetaData_t::AddMetaData(std::string(opt),
                          std::string(val),
                          std::string(err_msg)
                          );
  return def_val;
}

void commandline_option_end(int argc, char** argv){
  char help[]="--help";
  for(int i=0;i<argc;i++){
    if(!strcmp(argv[i],help)){
      std::cout<<"\n";
      exit(0);
    }
  }
}

void parse_command_line_options(int argc, char** argv) {
  // Read command line options.
  commandline_option_start(argc, argv);
  int omp =
      atoi(
          commandline_option(argc, argv, "-omp", "1", false,
                             "-omp  <int> = (1)    : Number of OpenMP threads"));

  int np;
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  // default value is number of processes

  std::stringstream temp_str;
  temp_str<<np;
  std::string str = temp_str.str();

  size_t N =
      (size_t)strtod(
          commandline_option(argc, argv, "-N", str.c_str(), false,
                             "-N    <int>          : Number of point sources"),
          NULL);

  size_t M =
      (size_t)strtod(
          commandline_option(argc, argv, "-M", "1", false,
                             "-M    <int>          : Number of points per octant"),
          NULL);

  int q =
      strtoul(
          commandline_option(argc, argv, "-q", "14", false,
                             "-q    <int> = (14)   : Chebyshev order (+ve integer)"),
          NULL,10);

  int d =
      strtoul(
          commandline_option(argc, argv, "-d", "15", false,
                             "-d    <int> = (15)   : Maximum tree depth"),
          NULL,10);

  double tol =
      strtod(
          commandline_option(argc, argv, "-tol", "1e-5", false,
                             "-tol <real> = (1e-5) : Tolerance for adaptive refinement")
          ,NULL);
  bool adap =
      (commandline_option(argc, argv, "-adap", NULL, false,
                          "-adap                : Adaptive tree refinement" )!=NULL);

  int tn =
      strtoul(
          commandline_option(argc, argv, "-tn", "1", false,
                             "-tn   <int> = (1)    : Number of time steps"),
          NULL,10);

  double dt =
      strtod(
          commandline_option(argc, argv,  "-dt",  "1e-2", false,
                             "-dt <real> = (0.1e-2) : Temporal resolution" ), NULL);

  bool cubic =
      (commandline_option(argc, argv, "-cubic", NULL, false,
                          "-cubic               : Cubic Interpolation used to evaluate tree values")!=NULL);

  int cuf =
      strtoul(
          commandline_option(argc, argv, "-cuf", "4", false,
                             "-cuf   <int> = (4)    : Upsampling factor used for cubic interpolation"),
          NULL,10);

  bool profile =
      (commandline_option(argc, argv, "-p", NULL, false,
                          "-p                  : Deactivate profiling")==NULL);

  int vtk_save_rate =
      atoi(
          commandline_option(argc, argv, "-vsr", "1", false,
                             "-vsr  <int> = (1)    : VTK files saving rate"));


  commandline_option_end(argc, argv);
  // =========================================================================
  // SIMULATION PARAMETERS
  // =========================================================================
  tbslas::SimConfig* sim_config       = tbslas::SimConfigSingleton::Instance();

  sim_config->total_num_timestep      = tn;
  sim_config->dt                      = dt;
  sim_config->vtk_order               = q;
  sim_config->vtk_save_rate           = vtk_save_rate;
  sim_config->use_cubic               = cubic;
  sim_config->cubic_upsampling_factor = cubic?cuf:0;
  sim_config->num_omp_threads         = omp;
  omp_set_num_threads(omp);
  // *************************************************************************
  // chebyshev tree
  // *************************************************************************
  sim_config->tree_num_point_sources      = N;
  sim_config->tree_num_points_per_octanct = M;
  sim_config->tree_chebyshev_order        = q;
  sim_config->tree_max_depth              = d;
  sim_config->tree_tolerance              = tol;
  sim_config->tree_adap                   = adap;
  // *************************************************************************
  // MISC
  // *************************************************************************
  sim_config->profile = profile;
}
