/* file: Gas.cpp

	This is the implementation file for the gas class, it will eventually be the object for
	each gas that populates a current gas model.

    To remind myself of the physics, I'm trying to compute stuff based
	on basic principles at start. 
	
	compile: g++ -Wall Gas.cpp -o gas
	
	Megan Sosey
	December 2011
*/

#include <iostream>
#include <cassert>
#include <cmath>
#include "Gas.h"

using std::cout;
using std::endl;

	
/*default constructor which makes a new Gas object and initializes the state and name */
Gas::Gas() {
			
	stateName=ATOMIC;
	gasName=H; //default to the simplest gas
	
}

//Alternate constructor which lets you specify everything
Gas::Gas(STATE state, GASNAME gas){

	stateName=state;
	gasName=gas;
}

Gas::~Gas(){}

int Gas::get_atomic_number(void){
	//return the atomic number for the gas, simplistic
	return int(gasName) + 1;
}

void Gas::setState(STATE state) {
	//set the state of the gas
	stateName=state;
}

STATE Gas::getState(void)  {
	//return the state of the gas
	return(stateName);
}


GASNAME Gas::getName(void)  {
	//return the gas name
	return(gasName);
}
	
float Gas::compute_orbital_radius(int n) {

	//return the radius for the energy level provided
	return (n*n*bohr_radius);
}

float Gas::compute_ionization_energy(int n){

	//compute the energy necessary to ionize hydrogen atom from the given level
	//E = R/n**2 where n = level   computed in ergs
	//bohr estimate => E= R Z**2 / n**2  where Z= the atomic number
	int atomic_no=get_atomic_number();
	return rydberg_R  * (atomic_no * atomic_no) / (n*n);	

}


float Gas::compute_ionization_wavelength(float energy) {

	//return the wavelength emmitted when an energy is released
	//according to E = hc/lambda
	return (plank_h * speed_c) / energy;
}

/* float Gas::compute_larmor_power(int Z){
	//compute the em radiation power(P) for an electrical charge q for electrostatic bremsstrahlung
	// P = 2 q**2 v**2 / 3 c**3
	// v = -Ze**2/ m(electron)l**2
	float Z = float(get_atomic_number());
	return (  Z* Z * e_charge * e_charge * e_charge * e_charge ) / (e_mass* e_mass *3. * speed_c*speed_c*speed_c);
	
}
 */
 
float Gas::compute_bohr_radius(void){

	//calculate the reduced Bohr radius for the atom?
	//compton wave is h/2pi * m * c
	//float compton_e = plank_h/ (2.* pi * e_mass * speed_c);
	//float compton_p = plank_h/ (2.* pi * p_mass * speed_c);
	//return  (compton_p + compton_e) / (2 * pi * fine_structure_a); //reduced bohr radius
	
    // the regular radius =  plank_h / (2.*pi * e_mass * speed_c * fine_structure_a);
	return plank_h / (2.*pi * e_mass * speed_c * fine_structure_a);
}

float Gas::compute_avg_photon_energy(float tstar) {
	// calculation the energy of the average photon at tstar temp using wiens law
	// wave(max) =0.29/T [cm] and energy is hc/wave
	float wave = 0.29 / tstar;
	return ( plank_h * speed_c )/wave; 
}

float Gas::compute_absorption_csection(void){

	//calculate the absorption cross section, approximate as the (reduced Bohr radius for the atom) ** 2
	//is this ok?
	float bohr = compute_bohr_radius();
	float atomic_no=float(get_atomic_number());
	float r = (bohr * atomic_no);
	return  (r * r);
}

double Gas::compute_stromgren_radius_thickness(double density){

	//calculate the thickness of the partially ionized skin
	double csection = compute_absorption_csection();
	return  1./(density * csection); //[cm]
	
}

float Gas::compute_h_recomb_rate(float tstar, float e_density, float p_density){

	//compute the Hydrogen recombination rate to form neutral hydrogen
	//dn/dt = -2e-13 * e_density * p_density (tstar / 1e4K) ^ (-3/4) per cm**3 per s
	float first = 2e-13 * e_density * p_density;
	float temp = pow(tstar / 1e4, 4./3.);
	return first * temp;
	
	
}

double Gas::compute_stromgren_equil(double n_density, float e_density, float tstar){
	
	//calculate the radius of the fully ionized stromgren sphere
	// Rs = (( 3. * n_density) / (4pi * hydrogen recomb rate * e_density**2) ) **1/3
	//n_density is the number of lyman continuum photons
	//e_density is the number of electrons
	
    double numer = (3. *  n_density);
	double h_recombination_rate = double(compute_h_recomb_rate(tstar, e_density, e_density));
	double denom = (4. * pi * h_recombination_rate * e_density * e_density);
	
	return pow(numer/denom, 1./3.) ; //[cm]
}

float Gas::compute_recombination_rate(float volume, float e_density, float p_density, float tstar){
	//compute the volume dependent recombination rate, per time (cm**-3 s**-1)
	// rate(volume) = hydrogen recombination rate * n(electron density) * n(proton density)
	double h_recombination_rate = double(compute_h_recomb_rate(tstar, e_density, p_density));
	return h_recombination_rate * e_density * p_density * volume;
}

float Gas::compute_lyman_series(int start, int finish){
	//return the lyman series emission wavelengths for the atom
	// wave = hc / (Einit - Efinal)
	float hc = plank_h * speed_c;
	float Einit = compute_ionization_energy(start);
	float Efinal = compute_ionization_energy(finish);
	
	return (hc / (Efinal - Einit));
}


const char* Gas::printState()  {
	const char* stateNames[] = {"Atomic","Ionized","Molecular"}; //these should be in same order as STATE list, simplistic

	//print the nice name of the state
	return(stateNames[stateName]);
}

const char* Gas::printName() {
    const char*  gasNames[] = {"Hydrogen","Helium"};	//same order as the gas list, simplistic
	//print the nice name of the gas
	return(gasNames[gasName]);
}
	
	
/* test driver for the Gas object, easier to keep it with the code
   just call the compiled program and it will run
*/

int main (void) {

	float n0_cloud=1e3; //density of the HI cloud [atoms/cm3]
	unsigned int n_levels[]={1,2}; //different levels in the atom, ground transition from 2s -> 1s
	float collisional_recomb_rate=0;
	float electron_density = 1e3; 
	float proton_density = 1e3;
	float recomb_time = 0;
	double Nlyman = 6e49 ; // the number of Lyman continuum photons generated by star per second
	double thickness= 0.;
	float tstar = 44.5e3; //temperature of star in Kelvins
	
	cout <<"\nCreating default hydrogen atom:\n" << endl;
	Gas  hGas;
	cout << hGas.printName() << "\t" << hGas.printState() << "\t" << hGas.get_atomic_number() <<endl;
	
	//Now lets try computing it for all the levels in n
	cout << "Ionization energy [eV] and wavelength [A] for levels in  " << hGas.printName() <<endl;
	for (unsigned int i=0;i< sizeof(n_levels)/sizeof(unsigned int); i++){
		float energy = hGas.compute_ionization_energy(n_levels[i]);
		float wave = hGas.compute_ionization_wavelength(energy) / angstroms;
		cout << n_levels[i] << "\t" << energy/eV << "\t" << wave <<endl;
	}

	//Now for some cross-sections and interaction	
	thickness = hGas.compute_stromgren_radius_thickness(n0_cloud) ;
	cout << "\nUsing a density of "<< n0_cloud << " atoms/cm**3, the Stromgren radius thickness is: " << thickness << "[cm]" <<endl;
	cout << "The absorption cross-section for the " << hGas.printName() << " atom was " << hGas.compute_absorption_csection() <<endl;
	
	//compute energy and radius the other way, using the temperature of the star
	cout << "Ionization energy [eV] and wavelength [A] for levels in  " << hGas.printName() <<endl;
		float energy = hGas.compute_avg_photon_energy(tstar);
		float wave = hGas.compute_ionization_wavelength(energy) / angstroms;
		cout <<energy/eV << "\t" << wave <<endl;
	
	cout << "\nStromgren radius in equilibrium: " << hGas.compute_stromgren_equil(Nlyman, electron_density, tstar) << "[cm]" <<endl;
	
	collisional_recomb_rate = hGas.compute_recombination_rate(n0_cloud, electron_density, proton_density, tstar);
	cout << "\nThe collisional recombination rate =  " << collisional_recomb_rate << " [per cm**-3 per s]" << endl;
	
	recomb_time =  electron_density / collisional_recomb_rate;
	cout << "\nThe recombination time is: "<< recomb_time << " seconds (" << recomb_time / 3600. << " hours)" <<endl;
	
	
	cout << "\nLyman alpha for n=2 " << hGas.compute_lyman_series(2,1) / angstroms << "[A]" <<endl;
	cout << "\nLyman beta for n=3 " << hGas.compute_lyman_series(3,1) / angstroms << "[A]" <<endl;
	
	cout <<"\n\nCreating  helium:\n";
	Gas  heGas(ATOMIC,HE);
	cout << heGas.printName() << "\t" << heGas.printState()<< "\t" << heGas.get_atomic_number() <<endl;

	//Now lets try computing it for all the levels in n
	cout << "Ionization energy [eV] and wavelength [A] for levels in  " << heGas.printName() <<endl;
	for (unsigned int i=0;i<sizeof(n_levels)/sizeof(unsigned int);i++){
		float energy = heGas.compute_ionization_energy(n_levels[i]);
		float wave = heGas.compute_ionization_wavelength(energy) / angstroms;
		cout << n_levels[i] << "\t" << energy/eV << "\t" << wave <<endl;
	}
	
	
	return(0);
	
}
