#ifndef	_PAW_BASIS_H_
#define _PAW_BASIS_H_
/*
   $Id$
*/


extern void paw_init_paw_basis(
        char*,
        char*,
        int ,
        int*,
        int*,
        float*
    );

extern void paw_end_paw_basis();

extern void  paw_solve_pseudo_orbitals();

extern void paw_generate_pseudo_orbitals();

extern void    paw_generate_projectors();

extern void    paw_generate_projectors_vanderbilt();

extern void    paw_generate_projectors_blochl();

extern void paw_print_basis_to_file(char* );

extern void paw_get_original_density(float *rho,float *rho_ps);

extern float* paw_get_pointer_paw_ps_density();

extern float* paw_get_pointer_paw_density();

extern float paw_get_paw_kinetic_energy();
extern float paw_get_r_orbital();

extern int paw_get_nbasis();
extern float paw_get_Zvalence();

extern int* paw_get_pointer_paw_l_array();

extern int* paw_get_pointer_paw_n_array();

extern int* paw_get_pointer_paw_n_ps_array();

extern float* paw_get_pointer_paw_e_array();

extern float** paw_get_pointer_paw_psi_array();

extern float** paw_get_pointer_paw_psi_prime_array();

extern float** paw_get_pointer_paw_psi_ps_array();

extern float** paw_get_pointer_paw_psi_ps_prime_array();

extern float** paw_get_pointer_paw_prj_ps_array();

extern float** paw_get_pointer_paw_prj_ps0_array();

extern int paw_get_max_i_r_orbital();

extern int paw_projectors_are_done();

extern float** paw_get_pointer_paw_psi_ps_unscr_array();

extern void paw_print_basis_information(FILE *fp);

extern void paw_print_basis_test_to_file(char* atom_name);
extern void paw_scattering_test(float e1,float e2,int number_points ,int l, float r );
#endif

