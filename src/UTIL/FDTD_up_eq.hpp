#ifndef PARALLEL_FDTD_UPEQ
#define PARALLEL_FDTD_UPEQ

#include <PML/parallelPML.hpp>
#include <UTIL/ml_consts.hpp>

namespace FDTDCompUpdateFxnReal
{
    typedef real_pgrid_ptr pgrid_ptr;
    /**
     * @brief      Updates the fields by taking the curl and assuming the k field is not there
     *
     * @param[in]  axList      parameters for the blas functions
     * @param[in]  axParams    temp array to store the axpy parameters
     * @param[in]  prefactors  temp array to store the prefactors
     * @param[in]  grid_i      shared_ptr to the field to be updated
     * @param[in]  grid_j      shared_ptr to the jth field with ijk notation: (i.e. if grid_i is Ey, grid_j is Hz)
     * @param[in]  grid_k      shared_ptr to the kth field with ijk notation: (i.e. if grid_i is Ex, grid_j is Hz)
     */
    void OneCompCurlJ (std::array<int,8>& axParams, std::array<double,2>& prefactors, pgrid_ptr grid_i, pgrid_ptr grid_j, pgrid_ptr grid_k);

    /**
     * @brief      Updates the fields by taking the curl and assuming the j field is not there
     *
     * @param[in]  axList      parameters for the blas functions
     * @param[in]  axParams    temp array to store the axpy parameters
     * @param[in]  prefactors  temp array to store the prefactors
     * @param[in]  grid_i      shared_ptr to the field to be updated
     * @param[in]  grid_j      shared_ptr to the jth field with ijk notation: (i.e. if grid_i is Ey, grid_j is Hz)
     * @param[in]  grid_k      shared_ptr to the kth field with ijk notation: (i.e. if grid_i is Ex, grid_j is Hz)
     */
    void OneCompCurlK (std::array<int,8>& axParams, std::array<double,2>& prefactors, pgrid_ptr grid_i, pgrid_ptr grid_j, pgrid_ptr grid_k);

    /**
     * @brief      Updates the fields by taking the curl and assuming that both fields are there
     *
     * @param[in]  axList      parameters for the blas functions
     * @param[in]  axParams    temp array to store the axpy parameters
     * @param[in]  prefactors  temp array to store the prefactors
     * @param[in]  grid_i      shared_ptr to the field to be updated
     * @param[in]  grid_j      shared_ptr to the jth field with ijk notation: (i.e. if grid_i is Ey, grid_j is Hz)
     * @param[in]  grid_k      shared_ptr to the kth field with ijk notation: (i.e. if grid_i is Ex, grid_j is Hz)
     */
    void TwoCompCurl (std::array<int,8>& axParams, std::array<double,2>& prefactors, pgrid_ptr grid_i, pgrid_ptr grid_j, pgrid_ptr grid_k);

    /**
     * @brief      Updates the Lorentzian Polarization Fields
     *
     * @param[in]  axList      List of update parameters for ?axpy_
     * @param[in]  axParams    temp array for ?axpy_ updates
     * @param[in]  grid_i      Electric field corresponding to polarization that will be updated
     * @param[in]  lorPi       vector of polarization fields (one for every mode)
     * @param[in]  prevLorPi   vector of polarization fields at the previous time step
     * @param[in]  jstore      pointer to scratch space to store the J values
     * @param[in]  obj         shared_ptr to the object that the field is interacting with
     */
    void UpdateLorPol(std::array<int,8>& axParams, pgrid_ptr grid_i, std::vector<pgrid_ptr> & lorPi, std::vector<pgrid_ptr> & prevLorPi, double* jstore, std::shared_ptr<Obj> obj);

    /**
     * @brief      Updates the Lorentzian Polarization Fields
     *
     * @param[in]  axList      List of update parameters for ?axpy_
     * @param[in]  axParams    temp array for ?axpy_ updates
     * @param[in]  grid_i      Magnetic field corresponding to polarization that will be updated
     * @param[in]  lorMi       vector of magnetization fields (one for every mode)
     * @param[in]  prevLorMi   vector of magnetization fields at the previous time step
     * @param[in]  jstore      pointer to scratch space to store the J values
     * @param[in]  obj         shared_ptr to the object that the field is interacting with
     */
    void UpdateLorMag(std::array<int,8>& axParams, pgrid_ptr grid_i, std::vector<pgrid_ptr> & lorMi, std::vector<pgrid_ptr> & prevLorMi, double* jstore, std::shared_ptr<Obj> obj);

    /**
     * @brief      Takes D field and moves it to the E field
     *
     * @param[in]  axList      List of update parameters for ?axpy_
     * @param[in]  axParams    temp array for ?axpy_ updates
     * @param[in]  grid_i      shared_ptr to the electric field corresponding to polarization that will be updated
     * @param[in]  lorPi       vector of shared_ptrs to the polarization fields (one for every mode)
     * @param[in]  prevLorPi   vector of shared_ptrs to the polarization fields at the previous time step
     * @param[in]  jstore      dummy vector for storing current polarization field
     * @param[in]  obj         shared_ptr to the object that the field is interacting with
     */
    void DtoE(std::array<int,8>& axParams, pgrid_ptr Di, pgrid_ptr Ei, std::vector<pgrid_ptr> & lorPi, std::shared_ptr<Obj> obj);


    /**
     * @brief      Takes D field and moves it to the E field
     *
     * @param[in]  axList      List of update parameters for ?axpy_
     * @param[in]  axParams    temp array for ?axpy_ updates
     * @param[in]  grid_i      shared_ptr to the electric field corresponding to polarization that will be updated
     * @param[in]  lorPi       vector of shared_ptrs to the polarization fields (one for every mode)
     * @param[in]  prevLorPi   vector of shared_ptrs to the polarization fields at the previous time step
     * @param[in]  jstore      dummy vector for storing current polarization field
     * @param[in]  obj         shared_ptr to the object that the field is interacting with
     */
    void BtoH(std::array<int,8>& axParams, pgrid_ptr Bi, pgrid_ptr Hi, std::vector<pgrid_ptr> & lorMi, std::shared_ptr<Obj> obj);

    /**
     * @brief      applies periodic boundary conditions to the field
     *
     * @param[in]  fUp      shared_ptr to the field that is being updated
     * @param[in]  k_point  vector representing the k_point of light
     * @param[in]  nx       right boundary of the field
     * @param[in]  ny       top boundary of the cell
     * @param[in]  dx       grid spacing in the x direction
     * @param[in]  dy       grid spacing in the y direction
     */
    void applyPBC(pgrid_ptr fUp, std::array<double,3> & k_point, int nx, int ny, int nz, int xmax, int ymax, int zmin, int zmax, double & dx, double & dy, double & dz);

    /**
     * @brief      applies periodic boundary conditions to the field
     *
     * @param[in]  fUp      shared_ptr to the field that is being updated
     * @param[in]  k_point  vector representing the k_point of light
     * @param[in]  nx       right boundary of the field
     * @param[in]  ny       top boundary of the cell
     * @param[in]  dx       grid spacing in the x direction
     * @param[in]  dy       grid spacing in the y direction
     */
    void applyPBC1Proc(pgrid_ptr fUp, std::array<double,3> & k_point, int nx, int ny, int nz, int xmax, int ymax, int zmin, int zmax, double & dx, double & dy, double & dz);

}

namespace FDTDCompUpdateFxnCplx
{
    typedef cplx_pgrid_ptr pgrid_ptr;

    /**
     * @brief      Updates the fields by taking the curl and assuming the k field is not there
     *
     * @param[in]  axList      parameters for the blas functions
     * @param[in]  axParams    temp array to store the axpy parameters
     * @param[in]  prefactors  temp array to store the prefactors
     * @param[in]  grid_i      shared_ptr to the field to be updated
     * @param[in]  grid_j      shared_ptr to the jth field with ijk notation: (i.e. if grid_i is Ey, grid_j is Hz)
     * @param[in]  grid_k      shared_ptr to the kth field with ijk notation: (i.e. if grid_i is Ex, grid_j is Hz)
     */
    void OneCompCurlJ (std::array<int,8>& axParams, std::array<double,2>& prefactors, pgrid_ptr grid_i, pgrid_ptr grid_j, pgrid_ptr grid_k);

    /**
     * @brief      Updates the fields by taking the curl and assuming the j field is not there
     *
     * @param[in]  axList      parameters for the blas functions
     * @param[in]  axParams    temp array to store the axpy parameters
     * @param[in]  prefactors  temp array to store the prefactors
     * @param[in]  grid_i      shared_ptr to the field to be updated
     * @param[in]  grid_j      shared_ptr to the jth field with ijk notation: (i.e. if grid_i is Ey, grid_j is Hz)
     * @param[in]  grid_k      shared_ptr to the kth field with ijk notation: (i.e. if grid_i is Ex, grid_j is Hz)
     */
    void OneCompCurlK (std::array<int,8>& axParams, std::array<double,2>& prefactors, pgrid_ptr grid_i, pgrid_ptr grid_j, pgrid_ptr grid_k);

    /**
     * @brief      Updates the fields by taking the curl and assuming that both fields are there
     *
     * @param[in]  axList      parameters for the blas functions
     * @param[in]  axParams    temp array to store the axpy parameters
     * @param[in]  prefactors  temp array to store the prefactors
     * @param[in]  grid_i      shared_ptr to the field to be updated
     * @param[in]  grid_j      shared_ptr to the jth field with ijk notation: (i.e. if grid_i is Ey, grid_j is Hz)
     * @param[in]  grid_k      shared_ptr to the kth field with ijk notation: (i.e. if grid_i is Ex, grid_j is Hz)
     */
    void TwoCompCurl (std::array<int,8>& axParams, std::array<double,2>& prefactors, pgrid_ptr grid_i, pgrid_ptr grid_j, pgrid_ptr grid_k);

    /**
     * @brief      Updates the Lorentzian Polarization Fields
     *
     * @param[in]  axList      List of update parameters for ?axpy_
     * @param[in]  axParams    temp array for ?axpy_ updates
     * @param[in]  grid_i      shared_ptr to the electric field corresponding to polarization that will be updated
     * @param[in]  lorPi       vector of shared_ptrs to the polarization fields (one for every mode)
     * @param[in]  prevLorPi   vector of shared_ptrs to the polarization fields at the previous time step
     * @param[in]  jstore      pointer to scratch space to store the J vectors
     * @param[in]  obj         shared_ptr to the object that the field is interacting with
     */
    void UpdateLorPol(std::array<int,8>& axParams, pgrid_ptr grid_i, std::vector<pgrid_ptr> & lorPi, std::vector<pgrid_ptr> & prevLorPi, cplx* jstore, std::shared_ptr<Obj> obj);

    /**
     * @brief      Updates the Lorentzian Polarization Fields
     *
     * @param[in]  axList      List of update parameters for ?axpy_
     * @param[in]  axParams    temp array for ?axpy_ updates
     * @param[in]  grid_i      Magnetic field corresponding to polarization that will be updated
     * @param[in]  lorMi       vector of magnetization fields (one for every mode)
     * @param[in]  prevLorMi   vector of magnetization fields at the previous time step
     * @param[in]  jstore      pointer to scratch space to store the J values
     * @param[in]  obj         shared_ptr to the object that the field is interacting with
     */
    void UpdateLorMag(std::array<int,8>& axParams, pgrid_ptr grid_i, std::vector<pgrid_ptr> & lorMi, std::vector<pgrid_ptr> & prevLorMi, cplx* jstore, std::shared_ptr<Obj> obj);

    /**
     * @brief      Takes D field and moves it to the E field
     *
     * @param[in]  axList      List of update parameters for ?axpy_
     * @param[in]  axParams    temp array for ?axpy_ updates
     * @param[in]  grid_i      shared_ptr to the electric field corresponding to polarization that will be updated
     * @param[in]  lorPi       vector of shared_ptrs to the polarization fields (one for every mode)
     * @param[in]  prevLorPi   vector of shared_ptrs to the polarization fields at the previous time step
     * @param[in]  jstore      dummy vector for storing current polarization field
     * @param[in]  obj         shared_ptr to the object that the field is interacting with
     */
    void DtoE(std::array<int,8>& axParams, pgrid_ptr Di, pgrid_ptr Ei, std::vector<pgrid_ptr> & lorPi, std::shared_ptr<Obj> obj);

    /**
     * @brief      Takes D field and moves it to the E field
     *
     * @param[in]  axList      List of update parameters for ?axpy_
     * @param[in]  axParams    temp array for ?axpy_ updates
     * @param[in]  grid_i      shared_ptr to the electric field corresponding to polarization that will be updated
     * @param[in]  lorPi       vector of shared_ptrs to the polarization fields (one for every mode)
     * @param[in]  prevLorPi   vector of shared_ptrs to the polarization fields at the previous time step
     * @param[in]  jstore      dummy vector for storing current polarization field
     * @param[in]  obj         shared_ptr to the object that the field is interacting with
     */
    void BtoH(std::array<int,8>& axParams, pgrid_ptr Bi, pgrid_ptr Hi, std::vector<pgrid_ptr> & lorMi, std::shared_ptr<Obj> obj);

    /**
     * @brief      applies periodic boundary conditions to the field
     *
     * @param[in]  fUp      shared_ptr to the field that is being updated
     * @param[in]  k_point  vector representing the k_point of light
     * @param[in]  nx       right boundary of the field
     * @param[in]  ny       top boundary of the cell
     * @param[in]  dx       grid spacing in the x direction
     * @param[in]  dy       grid spacing in the y direction
     */
    void applyPBC(pgrid_ptr fUp, std::array<double,3> & k_point, int nx, int ny, int nz, int xmax, int ymax, int zmin, int zmax, double & dx, double & dy, double & dz);

    /**
     * @brief      applies periodic boundary conditions to the field
     *
     * @param[in]  fUp      shared_ptr to the field that is being updated
     * @param[in]  k_point  vector representing the k_point of light
     * @param[in]  nx       right boundary of the field
     * @param[in]  ny       top boundary of the cell
     * @param[in]  dx       grid spacing in the x direction
     * @param[in]  dy       grid spacing in the y direction
     */
    void applyPBC1Proc(pgrid_ptr fUp, std::array<double,3> & k_point, int nx, int ny, int nz, int xmax, int ymax, int zmin, int zmax, double & dx, double & dy, double & dz);
}
#endif