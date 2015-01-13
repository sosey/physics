/* File: Gas.h

	This is the class declariation for the Gas implementation
	It holds all the information for each gas
	
	Megan Sosey
	December 2011
	
*/

#ifndef GAS_H
#define GAS_H
		
		
//Other useful constants for the calculations: 
const double plank_h = 6.6260755e-27; // [erg *s]
const double speed_c = 2.99792458e10; // [cm/s]
const double pi = 3.14159; //value for pi
const double fine_structure_a = 7.297353e-3; //constant
const double rydberg_R = 2.1798741e-11; //[ergs]
const double bohr_radius = 0.5291772e-10; //[cm] 			
const double eV = 1.602177e-12; //[ergs]
const double e_charge = 4.8032068e-10; //[esu]
const double joule = 1.e7; //[ergs]
const double e_mass = 9.1093897e-28;//[grams]
const double p_mass = 1.6726231e-24;//[grams]
const double angstroms = 1e-8; //[cm]
const double compton_p = 1.3214098446e-17; //[cm]
const double neutral_fraction = 0; //fully ionized n_e = n_p = n_H

//easier to read code lists 		
typedef enum STATE {ATOMIC, IONOIZED, MOLECULAR} ;
typedef enum GASNAME {H, HE} ; 

class Gas {
	
	public:
		//default constructor
		Gas();
		
		//alternate constructor
		Gas (STATE state, GASNAME gas);
		
		//destructor
		~Gas();
		
		//Associated Functions
		void setState(STATE stateName);		
		STATE getState(void);
		GASNAME getName(void);
		float compute_orbital_radius(int n);
		float compute_ionization_wavelength(float energy);
		float compute_ionization_energy(int n);
		//float compute_larmor_power(int Z);
		int get_atomic_number(void);
		double compute_stromgren_radius_thickness(double density);
		double compute_stromgren_equil(double n_density, float e_density, float tstar);
		float compute_absorption_csection(void);
		float compute_bohr_radius(void);
		float compute_recombination_rate(float volume, float e_density, float p_density, float tstar);
		float compute_lyman_series(int start, int finish);
		float compute_h_recomb_rate(float tstar, float e_density, float p_density);
		float compute_avg_photon_energy(float tstar);
		
		//pretty print information about the gas
  	    const char* printName (void);
		const char* printState(void);
		
	private:	
		STATE stateName;
		GASNAME gasName;	
		const char* stateNames[];
		const char* gasNames[];
};
	
#endif

