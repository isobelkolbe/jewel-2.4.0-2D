// Hydroinfo_MUSIC.cpp is a part of the MARTINI event generator.
// Copyright (C) 2009-2010 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains routines to read in hydro data from files and functions
// that return interpolated data at a given space-time point

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>

#include "data_struct.h"
#include "Hydroinfo_MUSIC.h"

using namespace std;


Hydroinfo_MUSIC::Hydroinfo_MUSIC() {
    boost_invariant = false;
}


Hydroinfo_MUSIC::~Hydroinfo_MUSIC() {
    if (boost_invariant) {
        lattice_2D.clear();
    } else {
        lattice_3D.clear();
        lattice_3D_ideal.clear();
    }
}


void Hydroinfo_MUSIC::readHydroData(int whichHydro, int nskip_tau_in) {
    // all hydro data is stored in tau steps (not t)
    // evolution is converted to tau when accessing the hydro data
    lattice_2D.clear();
    lattice_3D.clear();
    lattice_3D_ideal.clear();

    // read in setups of the hydro simulation
    hydroWhichHydro = whichHydro;
    if (hydroWhichHydro < 10) {
        ostringstream config_file;
        config_file << "results/music_input";
        ifstream configuration;
        configuration.open(config_file.str().c_str(), ios::in);
        if (!configuration) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << config_file.str() << endl;
            exit(1);
        }
        string temp1;
        string temp_name;
        while (!configuration.eof()) {
            getline(configuration, temp1);
            stringstream ss(temp1);
            ss >> temp_name;

            // read in grid information
            if (temp_name == "Initial_time_tau_0") {
                ss >> hydroTau0;
            } else if (temp_name == "Delta_Tau") {
                ss >> hydroDtau;
            } else if (temp_name == "X_grid_size_in_fm") {
                float temp;
                ss >> temp;
                hydroXmax = temp/2.;
            } else if (temp_name == "Grid_size_in_x") {
                ss >> ixmax;
            } else if (temp_name == "Eta_grid_size") {
                float temp;
                ss >> temp;
                hydro_eta_max = temp/2.;
            } else if (temp_name == "Grid_size_in_eta") {
                ss >> ietamax;
            } else if (temp_name == "output_evolution_every_N_timesteps") {
                ss >> nskip_tau;
            } else if (temp_name == "output_evolution_every_N_x") {
                ss >> nskip_x;
            } else if (temp_name == "output_evolution_every_N_eta") {
                ss >> nskip_eta;
            }
            // read in additioinal information
            if (temp_name == "Include_Rhob_Yes_1_No_0") {
                ss >> turn_on_rhob;
            } else if (temp_name == "Include_Shear_Visc_Yes_1_No_0") {
                ss >> turn_on_shear;
            } else if (temp_name == "Include_Bulk_Visc_Yes_1_No_0") {
                ss >> turn_on_bulk;
            } else if (temp_name == "turn_on_baryon_diffusion") {
                ss >> turn_on_diff;
            }
        }
        configuration.close();

        hydroDx = 2.*hydroXmax/(ixmax - 1.);
        hydroDeta = 2.*hydro_eta_max/(static_cast<float>(ietamax));
    }

    use_tau_eta_coordinate = 1;
    if (use_tau_eta_coordinate == 0) {
        cout << "Hydroinfo_MUSIC:: Warning hydro grid is set to "
             << "cartesian coordinates, please make sure this is correct!"
             << endl;
    }

    if (whichHydro == 6) {
        // 3+1D MUSIC hydro (Schenke, Jeon, Gale)
        cout << "Using 3+1D Jeon Schenke hydro reading data ..." << endl;
        boost_invariant = false;

        hydroDtau = hydroDtau*nskip_tau;
        hydroDx = 2.*hydroXmax/(ixmax - 1.)*nskip_x;
        hydroDeta = 2.*hydro_eta_max/(static_cast<float>(ietamax))*nskip_eta;

        ixmax = static_cast<int>(2.*hydroXmax/hydroDx + 0.001) + 1;
        ietamax = static_cast<int>(2.*hydro_eta_max/hydroDeta + 0.001);

        int n_eta = ietamax;
        int num_fluid_cell_trans = ixmax*ixmax;

        cout << "neta = " << n_eta
             << ", num_trans = " << ixmax << endl;

        // read in temperature, QGP fraction , flow velocity
        // The name of the evolution file: evolution_name
        string evolution_name = "evolution_xyeta.dat";
        string evolution_name_Wmunu =
            "results/evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
        string evolution_name_Pi = "results/evolution_bulk_pressure_xyeta.dat";
        cout << "Evolution file name = " << evolution_name << endl;

        std::FILE *fin;
        fin = std::fopen(evolution_name.c_str(), "rb");
        if (fin == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name << endl;
            exit(1);
        }

        std::FILE *fin1;
        fin1 = std::fopen(evolution_name_Wmunu.c_str(), "rb");
        if (fin1 == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name_Wmunu << endl;
            exit(1);
        }

        std::FILE *fin2;
        fin2 = std::fopen(evolution_name_Pi.c_str(), "rb");
        if (fin2 == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name_Pi << endl;
            exit(1);
        }

        fluidCell_3D newCell;
        int ik = 0;
        float T, QGPfrac, vx, vy, vz;
        float ux, uy, ueta;
        float pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33;
        float bulkPi, e_plus_P, cs2;
        int size = sizeof(float);
        while (true) {
            int status = 0;
            status = std::fread(&T, size, 1, fin);
            status += std::fread(&QGPfrac, size, 1, fin);
            status += std::fread(&vx, size, 1, fin);
            status += std::fread(&vy, size, 1, fin);
            status += std::fread(&vz, size, 1, fin);

            if (status != 5) {  // this is the end of file
                break;
            }

            int status_pi = 0;
            status_pi = std::fread(&pi00, size, 1, fin1);
            status_pi += std::fread(&pi01, size, 1, fin1);
            status_pi += std::fread(&pi02, size, 1, fin1);
            status_pi += std::fread(&pi03, size, 1, fin1);
            status_pi += std::fread(&pi11, size, 1, fin1);
            status_pi += std::fread(&pi12, size, 1, fin1);
            status_pi += std::fread(&pi13, size, 1, fin1);
            status_pi += std::fread(&pi22, size, 1, fin1);
            status_pi += std::fread(&pi23, size, 1, fin1);
            status_pi += std::fread(&pi33, size, 1, fin1);

            if (status_pi != 10) {
                cout << "Error:Hydroinfo_MUSIC::readHydroData: "
                     << "Wmunu file does not have the same number of "
                     << "fluid cells as the ideal file!" << endl;
                exit(1);
            }

            int status_bulkPi = 0;
            status_bulkPi = std::fread(&bulkPi, size, 1, fin2);
            status_bulkPi += std::fread(&e_plus_P, size, 1, fin2);
            status_bulkPi += std::fread(&cs2, size, 1, fin2);

            if (status_bulkPi != 3) {
                cout << "Error:Hydroinfo_MUSIC::readHydroData: "
                     << "bulkPi file does not have the same number of "
                     << "fluid cells as the ideal file!" << endl;
                exit(1);
            }

            int itau_idx = static_cast<int>(ik/(num_fluid_cell_trans*n_eta));
            ik++;

            int ieta = (static_cast<int>(ik/(num_fluid_cell_trans))
                        - n_eta*itau_idx);
            float eta = static_cast<float>(ieta)*hydroDeta - hydro_eta_max;

            // print out tau information
            float tau_local = hydroTau0 + itau_idx*hydroDtau;
            if ((ik-1)%(num_fluid_cell_trans*n_eta) == 0) {
                cout << "read in tau frame: " << itau_idx
                     << " tau_local = " << setprecision(3) << tau_local
                     << " fm ..."<< endl;
            }

            float v2 = vx*vx + vy*vy + vz*vz;
            if (v2 > 1.0) {
                if (T > 0.01) {
                    cerr << "[Hydroinfo_MUSIC::readHydroData:] Warning: "
                         << "v > 1! vx = " << vx << ", vy = " << vy
                         << ", vz = " << vz << ", T = " << T << endl;
                    exit(1);
                }
                ux = 0.0;
                uy = 0.0;
                ueta = 0.0;
            } else {
                float gamma = 1./sqrt(1. - v2);
                ux = gamma*vx;
                uy = gamma*vy;
                float uz = gamma*vz;
                ueta = -sinh(eta)*gamma + cosh(eta)*uz;
            }

            newCell.temperature = T;
            // convert vx and vy to longitudinal co-moving frame
            newCell.ux = ux;
            newCell.uy = uy;
            newCell.ueta = ueta;

            // pi^\mu\nu tensor
            newCell.pi00 = pi00;
            newCell.pi01 = pi01;
            newCell.pi02 = pi02;
            newCell.pi11 = pi11;
            newCell.pi12 = pi12;
            newCell.pi22 = pi22;
            newCell.pi33 = pi33;

            // bulk pressure
            if (T > 0.18) {
                // QGP phase prefactor is divided out here
                newCell.bulkPi = bulkPi/(15.*(1./3. - cs2)*e_plus_P);
            } else {
                newCell.bulkPi = bulkPi;   // [1/fm^4]
            }
            lattice_3D.push_back(newCell);
        }
        std::fclose(fin);
        std::fclose(fin1);
        std::fclose(fin2);
        cout << endl;
        cout << "number of fluid cells: " << lattice_3D.size() << endl;
    } else if (whichHydro == 8) {
        // event-by-event (2+1)-d MUSIC hydro from JF
        // there are two slices in medium in eta_s
        // one at eta_s = -15. and the other at eta_s = 0.0
        // only the medium at middle rapidity will be kept in the memory
        boost_invariant = true;
        cout << "Reading event-by-event hydro evolution data from JF ..."
             << endl;

        ixmax = static_cast<int>(2.*hydroXmax/hydroDx + 0.001);
        ietamax = 1;
        nskip_tau = nskip_tau_in;

        hydroDx *= nskip_x;
        hydroDtau *= nskip_tau;
        hydroDeta *= nskip_eta;

        int n_eta = 2;  // there are two slices in eta_s
        // number of fluid cell in the transverse plane
        int num_fluid_cell_trans = ixmax*ixmax;

        // read in hydro evolution
        string evolution_name = "results/evolution_xyeta.dat";
        string evolution_name_Wmunu =
                "results/evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
        string evolution_name_Pi = "results/evolution_bulk_pressure_xyeta.dat";

        std::FILE *fin;
        string evolution_file_name = evolution_name;
        cout << "Evolution file name = " << evolution_file_name << endl;
        fin = std::fopen(evolution_file_name.c_str(), "rb");

        if (fin == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_file_name << endl;
            exit(1);
        }

        std::FILE *fin1;
        fin1 = std::fopen(evolution_name_Wmunu.c_str(), "rb");
        if (fin1 == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name_Wmunu << endl;
            exit(1);
        }

        std::FILE *fin2;
        fin2 = std::fopen(evolution_name_Pi.c_str(), "rb");
        if (fin2 == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name_Pi << endl;
            exit(1);
        }

        int ik = 0;
        fluidCell_2D newCell;
        float T, QGPfrac, vx, vy, vz;
        float ux, uy, ueta;
        float pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33;
        float bulkPi, e_plus_P, cs2;
        int size = sizeof(float);
        while (true) {
            int status = 0;
            status = std::fread(&T, size, 1, fin);
            status += std::fread(&QGPfrac, size, 1, fin);
            status += std::fread(&vx, size, 1, fin);
            status += std::fread(&vy, size, 1, fin);
            status += std::fread(&vz, size, 1, fin);

            if (status != 5) {  // this is the end of file
                break;
            }

            int status_pi = 0;
            status_pi = std::fread(&pi00, size, 1, fin1);
            status_pi += std::fread(&pi01, size, 1, fin1);
            status_pi += std::fread(&pi02, size, 1, fin1);
            status_pi += std::fread(&pi03, size, 1, fin1);
            status_pi += std::fread(&pi11, size, 1, fin1);
            status_pi += std::fread(&pi12, size, 1, fin1);
            status_pi += std::fread(&pi13, size, 1, fin1);
            status_pi += std::fread(&pi22, size, 1, fin1);
            status_pi += std::fread(&pi23, size, 1, fin1);
            status_pi += std::fread(&pi33, size, 1, fin1);

            if (status_pi != 10) {
                cout << "Error:Hydroinfo_MUSIC::readHydroData: "
                     << "Wmunu file does not have the same number of "
                     << "fluid cells as the ideal file!" << endl;
                exit(1);
            }

            int status_bulkPi = 0;
            status_bulkPi = std::fread(&bulkPi, size, 1, fin2);
            status_bulkPi += std::fread(&e_plus_P, size, 1, fin2);
            status_bulkPi += std::fread(&cs2, size, 1, fin2);

            if (status_bulkPi != 3) {
                cout << "Error:Hydroinfo_MUSIC::readHydroData: "
                     << "bulkPi file does not have the same number of "
                     << "fluid cells as the ideal file!" << endl;
                exit(1);
            }

            int ieta_idx = static_cast<int>(ik/num_fluid_cell_trans) % n_eta;
            int itau_idx = static_cast<int>(ik/(num_fluid_cell_trans*n_eta));
            ik++;
            if (itau_idx%nskip_tau != 0)  // skip in tau
                continue;

            // print out tau information
            float tau_local = hydroTau0 + itau_idx*hydroDtau/nskip_tau;
            if ((ik-1)%(num_fluid_cell_trans*n_eta) == 0) {
                cout << "read in tau frame: " << itau_idx
                     << " tau_local = " << setprecision(3) << tau_local
                     << " fm ..."<< endl;
            }

            if (ieta_idx == (n_eta-1)) {
                // store the hydro medium at eta_s = 0.0
                float v2 = vx*vx + vy*vy + vz*vz;
                if (v2 > 1.0) {
                    cerr << "[Hydroinfo_MUSIC::readHydroData:] Error: "
                         << "v > 1! vx = " << vx << ", vy = " << vy
                         << ", vz = " << vz << ", T = " << T << endl;
                    if (T > 0.01) {
                        exit(1);
                    } else {
                        v2 = 0.0;
                    }
                }
                float gamma = 1./sqrt(1. - v2);
                ux = gamma*vx;
                uy = gamma*vy;
                ueta = gamma*vz;  // assuming eta = 0

                newCell.temperature = T;
                // convert vx and vy to longitudinal co-moving frame
                newCell.ux = ux;
                newCell.uy = uy;
                newCell.ueta = ueta;

                // pi^\mu\nu tensor
                newCell.pi00 = pi00;
                newCell.pi01 = pi01;
                newCell.pi02 = pi02;
                newCell.pi11 = pi11;
                newCell.pi12 = pi12;
                newCell.pi22 = pi22;
                newCell.pi33 = pi33;

                // bulk pressure
                if (T > 0.18) {
                    // QGP phase prefactor is divided out here
                    newCell.bulkPi = bulkPi/(15.*(1./3. - cs2)*e_plus_P);
                } else {
                    newCell.bulkPi = bulkPi;   // [1/fm^4]
                }
                lattice_2D.push_back(newCell);
            }
        }
        std::fclose(fin);
        std::fclose(fin1);
        std::fclose(fin2);
        cout << endl;
        cout << "number of fluid cells: " << lattice_2D.size() << endl;
    } else if (whichHydro == 9) {
        // event-by-event (2+1)-d MUSIC hydro
        // the output medium is at middle rapidity
        boost_invariant = true;
        cout << "Reading event-by-event hydro evolution data "
             << "from (2+1)D MUSIC ..." << endl;

        ietamax = 1;

        hydroDx *= nskip_x;
        hydroDtau *= nskip_tau;
        hydroDeta *= nskip_eta;

        nskip_tau = nskip_tau_in;
        ixmax = static_cast<int>(2.*hydroXmax/hydroDx + 0.001);

        int n_eta = 1;
        // number of fluid cell in the transverse plane
        int num_fluid_cell_trans = ixmax*ixmax;

        // read in hydro evolution
        string evolution_name = "results/evolution_xyeta.dat";
        string evolution_name_Wmunu =
                "results/evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
        string evolution_name_Pi = "results/evolution_bulk_pressure_xyeta.dat";

        std::FILE *fin;
        string evolution_file_name = evolution_name;
        cout << "Evolution file name = " << evolution_file_name << endl;
        fin = std::fopen(evolution_file_name.c_str(), "rb");

        if (fin == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_file_name << endl;
            exit(1);
        }

        std::FILE *fin1;
        fin1 = std::fopen(evolution_name_Wmunu.c_str(), "rb");
        if (fin1 == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name_Wmunu << endl;
            exit(1);
        }

        std::FILE *fin2;
        fin2 = std::fopen(evolution_name_Pi.c_str(), "rb");
        if (fin2 == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name_Pi << endl;
            exit(1);
        }

        int ik = 0;
        fluidCell_2D newCell;
        float T, QGPfrac, vx, vy, vz;
        float ux, uy, ueta;
        float pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33;
        float bulkPi, e_plus_P, cs2;
        int size = sizeof(float);
        while (true) {
            int status = 0;
            status = std::fread(&T, size, 1, fin);
            status += std::fread(&QGPfrac, size, 1, fin);
            status += std::fread(&vx, size, 1, fin);
            status += std::fread(&vy, size, 1, fin);
            status += std::fread(&vz, size, 1, fin);
            if (status != 5) {  // this is the end of file
                break;
            }

            float v2 = vx*vx + vy*vy + vz*vz;
            if (v2 > 1.) {
                cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                     << "v > 1! vx = " << vx << ", vy = " << vy
                     << ", vz = " << vz << endl;
                exit(1);
            }
            float gamma = 1./sqrt(1. - v2);
            ux = vx*gamma;
            uy = vy*gamma;
            ueta = vz*gamma;  // assuming at the eta = 0

            int status_pi = 0;
            status_pi = std::fread(&pi00, size, 1, fin1);
            status_pi += std::fread(&pi01, size, 1, fin1);
            status_pi += std::fread(&pi02, size, 1, fin1);
            status_pi += std::fread(&pi03, size, 1, fin1);
            status_pi += std::fread(&pi11, size, 1, fin1);
            status_pi += std::fread(&pi12, size, 1, fin1);
            status_pi += std::fread(&pi13, size, 1, fin1);
            status_pi += std::fread(&pi22, size, 1, fin1);
            status_pi += std::fread(&pi23, size, 1, fin1);
            status_pi += std::fread(&pi33, size, 1, fin1);

            if (status_pi != 10) {
                cout << "Error:Hydroinfo_MUSIC::readHydroData: "
                     << "Wmunu file does not have the same number of "
                     << "fluid cells as the ideal file!" << endl;
                exit(1);
            }

            int status_bulkPi = 0;
            status_bulkPi = std::fread(&bulkPi, size, 1, fin2);
            status_bulkPi += std::fread(&e_plus_P, size, 1, fin2);
            status_bulkPi += std::fread(&cs2, size, 1, fin2);

            if (status_bulkPi != 3) {
                cout << "Error:Hydroinfo_MUSIC::readHydroData: "
                     << "bulkPi file does not have the same number of "
                     << "fluid cells as the ideal file!" << endl;
                exit(1);
            }

            int ieta_idx = static_cast<int>(ik/num_fluid_cell_trans) % n_eta;
            int itau_idx = static_cast<int>(ik/(num_fluid_cell_trans*n_eta));
            ik++;
            if (itau_idx%nskip_tau != 0)  // skip in tau
                continue;

            // print out tau information
            float tau_local = hydroTau0 + itau_idx*hydroDtau/nskip_tau;
            if ((ik-1)%(num_fluid_cell_trans*n_eta) == 0) {
                cout << "read in tau frame: " << itau_idx
                     << " tau_local = " << setprecision(3) << tau_local
                     << " fm ..."<< endl;
            }

            if (ieta_idx == (n_eta-1)) {
                newCell.temperature = T;
                newCell.ux = ux;
                newCell.uy = uy;
                newCell.ueta = ueta;

                // pi^\mu\nu tensor
                newCell.pi00 = pi00;
                newCell.pi01 = pi01;
                newCell.pi02 = pi02;
                newCell.pi11 = pi11;
                newCell.pi12 = pi12;
                newCell.pi22 = pi22;
                newCell.pi33 = pi33;

                // bulk pressure
                if (T > 0.18) {
                    // QGP phase prefactor is divided out here
                    newCell.bulkPi = bulkPi/(15.*(1./3. - cs2)*e_plus_P);
                } else {
                    newCell.bulkPi = bulkPi;   // [1/fm^4]
                }
                lattice_2D.push_back(newCell);
            }
        }
        std::fclose(fin);
        std::fclose(fin1);
        std::fclose(fin2);
        cout << endl;
        cout << "number of fluid cells: " << lattice_2D.size() << endl;
    } else if (whichHydro == 10 || whichHydro == 11) {
        // new MUSIC hydro format (no grid)
        cout << "Using new MUSIC hydro format (no grid) reading data ..."
             << endl;
        if (whichHydro == 10) {
            boost_invariant = false;
        } else {
            boost_invariant = true;
        }

        // read in temperature and flow velocity
        // The name of the evolution file: evolution_name
        string evolution_name = "evolution_all_xyeta.dat";
        cout << "Evolution file name = " << evolution_name << endl;
        std::FILE *fin;
        fin = std::fopen(evolution_name.c_str(), "rb");
        if (fin == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name << endl;
            exit(1);
        }

        float header[16];
        int status = std::fread(&header, sizeof(float), 16, fin);
        if (status == 0) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Can not read the evolution file header" << endl;
            exit(1);
        }

        hydroTau0 = header[0];
        hydroDtau = header[1];
        ixmax = static_cast<int>(header[2]);
        hydroDx = header[3];
        hydroXmax = std::abs(header[4]);
        ietamax = static_cast<int>(header[8]);
        hydroDeta = header[9];
        hydro_eta_max = std::abs(header[10]);
        turn_on_rhob = static_cast<int>(header[11]);
        turn_on_shear = static_cast<int>(header[12]);
        turn_on_bulk = static_cast<int>(header[13]);
        turn_on_diff = static_cast<int>(header[14]);
        const int nVar_per_cell = static_cast<int>(header[15]);

        float cell_info[nVar_per_cell];

        int itau_max = 0;
        fluidCell_3D_ideal zeroCell;
        zeroCell.itau = 0;
        zeroCell.ix = 0;
        zeroCell.iy = 0;
        zeroCell.ieta = 0;
        zeroCell.temperature = 0.;
        zeroCell.ed = 0.;
        zeroCell.pressure = 0.;
        zeroCell.ux = 0.;
        zeroCell.uy = 0.;
        zeroCell.uz = 0.;
        lattice_3D_ideal.push_back(zeroCell);
        int ik = 0;
        while (true) {
            status = 0;
            status = std::fread(&cell_info, sizeof(float), nVar_per_cell, fin);
            if (status == 0) break;
            if (status != nVar_per_cell) {
                cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                     << "the evolution file format is not correct" << endl;
                exit(1);
            }

            if (itau_max < static_cast<int>(cell_info[0]))
                itau_max = static_cast<int>(cell_info[0]);
            fluidCell_3D_ideal newCell;
            newCell.itau = static_cast<int>(cell_info[0]);
            newCell.ix   = static_cast<int>(cell_info[1]);
            newCell.iy   = static_cast<int>(cell_info[2]);
            newCell.ieta = static_cast<int>(cell_info[3]);
            newCell.temperature = cell_info[6];
            newCell.ed = cell_info[4];
            newCell.pressure = cell_info[5];
            newCell.ux = cell_info[8];
            newCell.uy = cell_info[9];
            newCell.uz = cell_info[10];
            lattice_3D_ideal.push_back(newCell);
            ik++;
            if (ik%50000 == 0)
                cout << "o" << flush;
        }
        cout << endl;
        std::fclose(fin);
        itaumax = itau_max;
        // create the index map
        long long ncells = (itaumax + 1)*ixmax*ixmax*ietamax;
        idx_map_.resize(ncells, 0);
        for (int i = 0; i < lattice_3D_ideal.size(); i++) {
            const auto cell_i = lattice_3D_ideal[i];
            int cell_idx = (
                (  (cell_i.itau*ietamax + cell_i.ieta)*ixmax
                 + cell_i.iy)*ixmax + cell_i.ix);
            idx_map_[cell_idx] = i;
        }
        hydroTauMax = hydroTau0 + hydroDtau*itaumax;
    } else if (whichHydro == 13) {
        // new MUSIC hydro format (no grid)
        cout << "Using new MUSIC hydro format (no grid) "
             << "reading data with EM fields ..." << endl;
        boost_invariant = false;

        // read in temperature and flow velocity
        // The name of the evolution file: evolution_name
        string evolution_name = "results/evolution_all_xyeta_withEMfields.dat";
        cout << "Evolution file name = " << evolution_name << endl;
        std::FILE *fin;
        fin = std::fopen(evolution_name.c_str(), "rb");
        if (fin == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name << endl;
            exit(1);
        }

        float header[16];
        int status = std::fread(&header, sizeof(float), 16, fin);
        if (status == 0) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Can not read the evolution file header" << endl;
            exit(1);
        }

        hydroTau0 = header[0];
        hydroDtau = header[1];
        ixmax = static_cast<int>(header[2]);
        hydroDx = header[3];
        hydroXmax = std::abs(header[4]);
        ietamax = static_cast<int>(header[8]);
        hydroDeta = header[9];
        hydro_eta_max = std::abs(header[10]);
        turn_on_rhob = static_cast<int>(header[11]);
        turn_on_shear = static_cast<int>(header[12]);
        turn_on_bulk = static_cast<int>(header[13]);
        turn_on_diff = static_cast<int>(header[14]);
        const int nVar_per_cell = static_cast<int>(header[15]);

        float cell_info[nVar_per_cell];

        int itau_max = 0;
        fluidCell_3D_ideal_with_EM zeroCell;
        zeroCell.itau = 0;
        zeroCell.ix = 0;
        zeroCell.iy = 0;
        zeroCell.ieta = 0;
        zeroCell.temperature = 0.;
        zeroCell.ed = 0.;
        zeroCell.pressure = 0.;
        zeroCell.ux = 0.;
        zeroCell.uy = 0.;
        zeroCell.uz = 0.;
        zeroCell.Ex = 0.;
        zeroCell.Ey = 0.;
        zeroCell.Ez = 0.;
        zeroCell.Bx = 0.;
        zeroCell.By = 0.;
        zeroCell.Bz = 0.;
        lattice_3D_ideal_withEM_.push_back(zeroCell);
        int ik = 0;
        while (true) {
            status = 0;
            status = std::fread(&cell_info, sizeof(float), nVar_per_cell, fin);
            if (status == 0) break;
            if (status != nVar_per_cell) {
                cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                     << "the evolution file format is not correct" << endl;
                exit(1);
            }

            if (itau_max < static_cast<int>(cell_info[0]))
                itau_max = static_cast<int>(cell_info[0]);
            fluidCell_3D_ideal_with_EM newCell;
            newCell.itau = static_cast<int>(cell_info[0]);
            newCell.ix   = static_cast<int>(cell_info[1]);
            newCell.iy   = static_cast<int>(cell_info[2]);
            newCell.ieta = static_cast<int>(cell_info[3]);
            newCell.temperature = cell_info[6];
            newCell.ed = cell_info[4];
            newCell.pressure = cell_info[5];
            newCell.ux = cell_info[8];
            newCell.uy = cell_info[9];
            newCell.uz = cell_info[10];
            newCell.Ex = cell_info[11];
            newCell.Ey = cell_info[12];
            newCell.Ez = cell_info[13];
            newCell.Bx = cell_info[14];
            newCell.By = cell_info[15];
            newCell.Bz = cell_info[16];
            lattice_3D_ideal_withEM_.push_back(newCell);
            ik++;
            if (ik%50000 == 0)
                cout << "o" << flush;
        }
        cout << endl;
        std::fclose(fin);
        itaumax = itau_max;
        // create the index map
        long long ncells = (itaumax + 1)*ixmax*ixmax*ietamax;
        idx_map_.resize(ncells, 0);
        for (int i = 0; i < lattice_3D_ideal_withEM_.size(); i++) {
            const auto cell_i = lattice_3D_ideal_withEM_[i];
            int cell_idx = (
                (  (cell_i.itau*ietamax + cell_i.ieta)*ixmax
                 + cell_i.iy)*ixmax + cell_i.ix);
            idx_map_[cell_idx] = i;
        }
        hydroTauMax = hydroTau0 + hydroDtau*itaumax;
    } else {
        cout << "Hydroinfo_MUSIC:: This option is obsolete! whichHydro = "
             << whichHydro << endl;
        exit(1);
    }

    // One final step for easy automation of MARTINI:
    // hydroTauMax is reset for the case where writing to evolution.dat
    // ended early (due to all cells freezing out):
    if (whichHydro == 6) {
        hydroTauMax = (
            hydroTau0 + hydroDtau*static_cast<int>(
                        static_cast<float>(lattice_3D.size())
                        /((2.*hydroXmax/hydroDx+1.)*(2.*hydroXmax/hydroDx+1.)
                        *2.*(hydro_eta_max/hydroDeta))));
        itaumax = static_cast<int>((hydroTauMax-hydroTau0)/hydroDtau+0.001);
    }
    if (whichHydro == 8 || whichHydro == 9) {
        hydroTauMax = (
            hydroTau0 + hydroDtau*static_cast<int>(
                        static_cast<float>(lattice_2D.size())
                        /((2.*hydroXmax/hydroDx)*(2.*hydroXmax/hydroDx)) - 1));
        itaumax = static_cast<int>((hydroTauMax - hydroTau0)/hydroDtau);
    }

    cout << "hydro_tau0 = " << hydroTau0 << " fm"<< endl;
    cout << "hydro_tau_max = " << hydroTauMax << " fm" << endl;
    cout << "hydry_dtau = " << hydroDtau << " fm" << endl;
    cout << "hydro_Xmax = " << hydroXmax << " fm" << endl;
    cout << "hydro_dx = " << hydroDx << " fm" << endl;
    cout << "hydro_eta_max = " << hydro_eta_max << endl;
    cout << "hydro_deta = " << hydroDeta << endl;
}


void Hydroinfo_MUSIC::getHydroValues(float x, float y,
                                     float z, float t, fluidCell* info) {
// For interpolation of evolution files in tau-eta coordinates. Only the
// reading of MUSIC's evolution_xyeta.dat file is implemented here.
// For simplicity, hydro_eta_max refers to MUSIC's eta_size, and similarly for
// hydroDeta; however, x, y, z, and t are as usual to stay compatible with
// MARTINI.
    float tau, eta;
    if (use_tau_eta_coordinate == 1) {
        if (t*t > z*z) {
            tau = sqrt(t*t-z*z);
            eta = 0.5*log((t+z)/(t-z));
        } else {
            tau = 0.;
            eta = 0.;
        }
    } else {
        // if the medium is given in cartesian coordinates
        // set tau and eta to t and z
        tau = t;
        eta = z;
    }

    int itau = static_cast<int>((tau-hydroTau0)/hydroDtau + 0.0001);
    int ix   = static_cast<int>((hydroXmax+x)/hydroDx + 0.0001);
    int iy   = static_cast<int>((hydroXmax+y)/hydroDx + 0.0001);
    int ieta = static_cast<int>((hydro_eta_max+eta)/hydroDeta + 0.0001);

    float taufrac = (tau - hydroTau0)/hydroDtau - static_cast<float>(itau);
    float xfrac   = (x - (static_cast<float>(ix)*hydroDx - hydroXmax))/hydroDx;
    float yfrac   = (y - (static_cast<float>(iy)*hydroDx - hydroXmax))/hydroDx;
    float etafrac = (eta/hydroDeta - static_cast<float>(ieta)
                     + 0.5*static_cast<float>(ietamax));

    if (boost_invariant) {
        ieta = 0;
        etafrac = 0.;
    }

    if (ix < 0 || ix >= ixmax) {
        cout << "[Hydroinfo_MUSIC::getHydroValues]: "
             << "WARNING - x out of range x=" << x
             << ", ix=" << ix << ", ixmax=" << ixmax << endl;
        cout << "x=" << x << " y=" << y << " eta=" << eta
             << " ix=" << ix << " iy=" << iy << " ieta=" << ieta << endl;
        cout << "t=" << t << " tau=" << tau
             << " itau=" << itau << " itaumax=" << itaumax << endl;

        info->temperature = 0.0;
        info->vx = 0.0;
        info->vy = 0.0;
        info->vz = 0.0;
        return;
    }
    if (iy < 0 || iy >= ixmax) {
        cout << "[Hydroinfo_MUSIC::getHydroValues]: "
             << "WARNING - y out of range, y=" << y << ", iy="  << iy
             << ", iymax=" << ixmax << endl;
        cout << "x=" << x << " y=" << y << " eta=" << eta
             << " ix=" << ix << " iy=" << iy << " ieta=" << ieta << endl;
        cout << "t=" << t << " tau=" << tau
             << " itau=" << itau << " itaumax=" << itaumax << endl;

        info->temperature = 0.0;
        info->vx = 0.0;
        info->vy = 0.0;
        info->vz = 0.0;
        return;
    }
    if (itau < 0 || itau > itaumax) {
        cout << "[Hydroinfo_MUSIC::getHydroValues]: WARNING - "
             << "tau out of range, itau=" << itau << ", itaumax=" << itaumax
             << endl;
        cout << "[Hydroinfo_MUSIC::getHydroValues]: tau= " << tau
             << ", hydroTauMax = " << hydroTauMax
             << ", hydroDtau = " << hydroDtau << endl;

        info->temperature = 0.0;
        info->vx = 0.0;
        info->vy = 0.0;
        info->vz = 0.0;
        return;
    }
    if (ieta < 0 || ieta >= ietamax) {
        cout << "[Hydroinfo_MUSIC::getHydroValues]: WARNING - "
             << "eta out of range, ieta=" << ieta << ", ietamax=" << ietamax
             << endl;
        info->temperature = 0.0;
        info->vx = 0.0;
        info->vy = 0.0;
        info->vz = 0.0;
        return;
    }

    // The array of positions on the 4-dimensional rectangle:
    int position[2][2][2][2];
    for (int ipx = 0; ipx < 2; ipx++) {
        int px;
        if (ipx == 0 || ix == ixmax-1) {
            px = ix;
        } else {
            px = ix + 1;
        }
        for (int ipy = 0; ipy < 2; ipy++) {
            int py;
            if (ipy == 0 || iy == ixmax-1) {
                py = iy;
            } else {
                py = iy + 1;
            }
            for (int ipeta = 0; ipeta < 2; ipeta++) {
                int peta;
                if (ipeta == 0 || ieta == ietamax-1) {
                    peta = ieta;
                } else {
                    peta = ieta + 1;
                }
                for (int iptau = 0; iptau < 2; iptau++) {
                    int ptau;
                    if (iptau == 0 || itau == itaumax) {
                        ptau = itau;
                    } else {
                        ptau = itau + 1;
                    }
                    position[ipx][ipy][ipeta][iptau] = (
                                px + ixmax*(py + ixmax*(peta + ietamax*ptau)));
                }
            }
        }
    }

    // And now, the interpolation:
    float T = 0.0;
    float ed = 0.;
    float p = 0.;
    float vx = 0.0;
    float vy = 0.0;
    float vz = 0.0;
    float ux = 0.0;
    float uy = 0.0;
    float uz = 0.0;
    float ueta = 0.0;
    float pi00 = 0.0;
    float pi01 = 0.0;
    float pi02 = 0.0;
    float pi03 = 0.0;
    float pi11 = 0.0;
    float pi12 = 0.0;
    float pi13 = 0.0;
    float pi22 = 0.0;
    float pi23 = 0.0;
    float pi33 = 0.0;
    float bulkPi = 0.0;

    fluidCell_2D *HydroCell_2D_ptr1, *HydroCell_2D_ptr2;
    fluidCell_3D *HydroCell_3D_ptr1, *HydroCell_3D_ptr2;
    fluidCell_3D_ideal *HydroCell_3D_ideal_ptr1, *HydroCell_3D_ideal_ptr2;
    for (int iptau = 0; iptau < 2; iptau++) {
        float taufactor;
        if (iptau == 0)
            taufactor = 1. - taufrac;
        else
            taufactor = taufrac;
        for (int ipeta = 0; ipeta < 2; ipeta++) {
            float etafactor;
            if (ipeta == 0)
                etafactor = 1. - etafrac;
            else
                etafactor = etafrac;
            for (int ipy = 0; ipy < 2; ipy++) {
                float yfactor;
                if (ipy == 0)
                    yfactor = 1. - yfrac;
                else
                    yfactor = yfrac;

                float prefrac = yfactor*etafactor*taufactor;

                if (hydroWhichHydro == 8 || hydroWhichHydro == 9) {
                    HydroCell_2D_ptr1 = (
                            &lattice_2D[position[0][ipy][ipeta][iptau]]);
                    HydroCell_2D_ptr2 = (
                            &lattice_2D[position[1][ipy][ipeta][iptau]]);
                    T += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->temperature
                                  + xfrac*HydroCell_2D_ptr2->temperature);
                    ux += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->ux
                                    + xfrac*HydroCell_2D_ptr2->ux);
                    uy += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->uy
                                    + xfrac*HydroCell_2D_ptr2->uy);
                    ueta += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->ueta
                                    + xfrac*HydroCell_2D_ptr2->ueta);
                    if (turn_on_shear == 1) {
                        pi00 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi00
                                        + xfrac*HydroCell_2D_ptr2->pi00);
                        pi01 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi01
                                        + xfrac*HydroCell_2D_ptr2->pi01);
                        pi02 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi02
                                        + xfrac*HydroCell_2D_ptr2->pi02);
                        pi11 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi11
                                        + xfrac*HydroCell_2D_ptr2->pi11);
                        pi12 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi12
                                        + xfrac*HydroCell_2D_ptr2->pi12);
                        pi22 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi22
                                        + xfrac*HydroCell_2D_ptr2->pi22);
                        pi33 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi33
                                        + xfrac*HydroCell_2D_ptr2->pi33);
                    }
                    if (turn_on_bulk == 1) {
                        bulkPi += (
                            prefrac*((1. - xfrac)*HydroCell_2D_ptr1->bulkPi
                            + xfrac*HydroCell_2D_ptr2->bulkPi));
                    }
                } else if (hydroWhichHydro == 6) {
                    HydroCell_3D_ptr1 = (
                            &lattice_3D[position[0][ipy][ipeta][iptau]]);
                    HydroCell_3D_ptr2 = (
                            &lattice_3D[position[1][ipy][ipeta][iptau]]);
                    T += prefrac*(
                            (1. - xfrac)*HydroCell_3D_ptr1->temperature
                            + xfrac*HydroCell_3D_ptr2->temperature);
                    ux += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->ux
                                   + xfrac*HydroCell_3D_ptr2->ux);
                    uy += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->uy
                                   + xfrac*HydroCell_3D_ptr2->uy);
                    ueta += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->ueta
                                     + xfrac*HydroCell_3D_ptr2->ueta);
                    if (turn_on_shear == 1) {
                        pi00 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi00
                                         + xfrac*HydroCell_3D_ptr2->pi00);
                        pi01 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi01
                                         + xfrac*HydroCell_3D_ptr2->pi01);
                        pi02 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi02
                                         + xfrac*HydroCell_3D_ptr2->pi02);
                        pi03 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi03
                                         + xfrac*HydroCell_3D_ptr2->pi03);
                        pi11 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi11
                                         + xfrac*HydroCell_3D_ptr2->pi11);
                        pi12 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi12
                                         + xfrac*HydroCell_3D_ptr2->pi12);
                        pi13 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi13
                                         + xfrac*HydroCell_3D_ptr2->pi13);
                        pi22 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi22
                                         + xfrac*HydroCell_3D_ptr2->pi22);
                        pi23 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi23
                                         + xfrac*HydroCell_3D_ptr2->pi23);
                        pi33 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi33
                                         + xfrac*HydroCell_3D_ptr2->pi33);
                    }
                    if (turn_on_bulk == 1) {
                        bulkPi += (
                            prefrac*((1. - xfrac)*HydroCell_3D_ptr1->bulkPi
                            + xfrac*HydroCell_3D_ptr2->bulkPi));
                    }
                } else if (hydroWhichHydro == 10 || hydroWhichHydro == 11) {
                    HydroCell_3D_ideal_ptr1 = (
                        &lattice_3D_ideal[idx_map_[position[0][ipy][ipeta][iptau]]]);
                    HydroCell_3D_ideal_ptr2 = (
                        &lattice_3D_ideal[idx_map_[position[1][ipy][ipeta][iptau]]]);
                    T += prefrac*(
                        (1. - xfrac)*HydroCell_3D_ideal_ptr1->temperature
                        + xfrac*HydroCell_3D_ideal_ptr2->temperature);
                    ed += prefrac*(
                        (1. - xfrac)*HydroCell_3D_ideal_ptr1->ed
                        + xfrac*HydroCell_3D_ideal_ptr2->ed);
                    p += prefrac*(
                        (1. - xfrac)*HydroCell_3D_ideal_ptr1->pressure
                        + xfrac*HydroCell_3D_ideal_ptr2->pressure);
                    ux += prefrac*((1. - xfrac)*HydroCell_3D_ideal_ptr1->ux
                                    + xfrac*HydroCell_3D_ideal_ptr2->ux);
                    uy += prefrac*((1. - xfrac)*HydroCell_3D_ideal_ptr1->uy
                                    + xfrac*HydroCell_3D_ideal_ptr2->uy);
                    uz += prefrac*((1. - xfrac)*HydroCell_3D_ideal_ptr1->uz
                                   + xfrac*HydroCell_3D_ideal_ptr2->uz);
                }
            }
        }
    }

    if (hydroWhichHydro == 10) {
        float ut = sqrt(1. + ux*ux + uy*uy + uz*uz);
        vx = ux/ut;
        vy = uy/ut;
        vz = uz/ut;
    } else {
        float eta_local = 0.5*log((t + z)/(t - z));
        float sinh_eta = sinh(eta_local);
        float cosh_eta = cosh(eta_local);
        float utau = sqrt(1. + ux*ux + uy*uy + ueta*ueta);
        float uz = utau*sinh_eta + ueta*cosh_eta;
        float ut = utau*cosh_eta + ueta*sinh_eta;
        vx = ux/ut;
        vy = uy/ut;
        vz = uz/ut;
    }

    info->temperature = T;
    info->vx = vx;
    info->vy = vy;
    info->vz = vz;

    if (hydroWhichHydro < 10) {
        info->ed = 1.0;                 // pi's are already divided by e+P
        info->sd = 0.0;
        info->pressure = 0.0;
    } else {
        info->ed = ed;
        info->sd = (ed + p)/(T + 1e-16);
        info->pressure = p;
    }


    info->pi[0][0] = pi00;
    info->pi[0][1] = pi01;
    info->pi[0][2] = pi02;
    info->pi[0][3] = pi03;
    info->pi[1][0] = pi01;
    info->pi[1][1] = pi11;
    info->pi[1][2] = pi12;
    info->pi[1][3] = pi13;
    info->pi[2][0] = pi02;
    info->pi[2][1] = pi12;
    info->pi[2][2] = pi22;
    info->pi[2][3] = pi23;
    info->pi[3][0] = pi03;
    info->pi[3][1] = pi13;
    info->pi[3][2] = pi23;
    info->pi[3][3] = pi33;

    info->bulkPi = bulkPi;
    return;
}


void Hydroinfo_MUSIC::getHydroValuesWithEMFields(
        float x, float y, float z, float t, fluidCellIdealWithEM* info) {
    float tau, eta;
    if (t*t > z*z) {
        tau = sqrt(t*t - z*z);
        eta = 0.5*log((t + z)/(t - z));
    } else {
        tau = 0.;
        eta = 0.;
    }
    cout << "tau = " << tau << ", eta = " << eta << endl;

    int itau = static_cast<int>((tau-hydroTau0)/hydroDtau + 0.0001);
    int ix   = static_cast<int>((hydroXmax+x)/hydroDx + 0.0001);
    int iy   = static_cast<int>((hydroXmax+y)/hydroDx + 0.0001);
    int ieta = static_cast<int>((hydro_eta_max+eta)/hydroDeta + 0.0001);

    float taufrac = (tau - hydroTau0)/hydroDtau - static_cast<float>(itau);
    float xfrac   = (x - (static_cast<float>(ix)*hydroDx - hydroXmax))/hydroDx;
    float yfrac   = (y - (static_cast<float>(iy)*hydroDx - hydroXmax))/hydroDx;
    float etafrac = (eta/hydroDeta - static_cast<float>(ieta)
                     + 0.5*static_cast<float>(ietamax));

    if (boost_invariant) {
        ieta = 0;
        etafrac = 0.;
    }

    // initialize the info cell with 0
    info->temperature = 0.0;
    info->ed = 0.0;
    info->pressure = 0.0;
    info->sd = 0.0;
    info->vx = 0.0;
    info->vy = 0.0;
    info->vz = 0.0;
    info->Ex = 0.0;
    info->Ey = 0.0;
    info->Ez = 0.0;
    info->Bx = 0.0;
    info->By = 0.0;
    info->Bz = 0.0;

    if (ix < 0 || ix >= ixmax) {
        cout << "[Hydroinfo_MUSIC::getHydroValues]: "
             << "WARNING - x out of range x=" << x
             << ", ix=" << ix << ", ixmax=" << ixmax << endl;
        cout << "x=" << x << " y=" << y << " eta=" << eta
             << " ix=" << ix << " iy=" << iy << " ieta=" << ieta << endl;
        cout << "t=" << t << " tau=" << tau
             << " itau=" << itau << " itaumax=" << itaumax << endl;
        return;
    }
    if (iy < 0 || iy >= ixmax) {
        cout << "[Hydroinfo_MUSIC::getHydroValues]: "
             << "WARNING - y out of range, y=" << y << ", iy="  << iy
             << ", iymax=" << ixmax << endl;
        cout << "x=" << x << " y=" << y << " eta=" << eta
             << " ix=" << ix << " iy=" << iy << " ieta=" << ieta << endl;
        cout << "t=" << t << " tau=" << tau
             << " itau=" << itau << " itaumax=" << itaumax << endl;
        return;
    }
    if (itau < 0 || itau > itaumax) {
        cout << "[Hydroinfo_MUSIC::getHydroValues]: WARNING - "
             << "tau out of range, itau=" << itau << ", itaumax=" << itaumax
             << endl;
        cout << "[Hydroinfo_MUSIC::getHydroValues]: tau= " << tau
             << ", hydroTauMax = " << hydroTauMax
             << ", hydroDtau = " << hydroDtau << endl;
        return;
    }
    if (ieta < 0 || ieta >= ietamax) {
        cout << "[Hydroinfo_MUSIC::getHydroValues]: WARNING - "
             << "eta out of range, ieta=" << ieta << ", ietamax=" << ietamax
             << endl;
        return;
    }

    // The array of positions on the 4-dimensional rectangle:
    int position[2][2][2][2];
    for (int ipx = 0; ipx < 2; ipx++) {
        int px;
        if (ipx == 0 || ix == ixmax-1) {
            px = ix;
        } else {
            px = ix + 1;
        }
        for (int ipy = 0; ipy < 2; ipy++) {
            int py;
            if (ipy == 0 || iy == ixmax-1) {
                py = iy;
            } else {
                py = iy + 1;
            }
            for (int ipeta = 0; ipeta < 2; ipeta++) {
                int peta;
                if (ipeta == 0 || ieta == ietamax-1) {
                    peta = ieta;
                } else {
                    peta = ieta + 1;
                }
                for (int iptau = 0; iptau < 2; iptau++) {
                    int ptau;
                    if (iptau == 0 || itau == itaumax) {
                        ptau = itau;
                    } else {
                        ptau = itau + 1;
                    }
                    position[ipx][ipy][ipeta][iptau] = (
                                px + ixmax*(py + ixmax*(peta + ietamax*ptau)));
                }
            }
        }
    }

    // And now, the interpolation:
    float T = 0.0;
    float ed = 0.;
    float p = 0.;
    float vx = 0.0;
    float vy = 0.0;
    float vz = 0.0;
    float ux = 0.0;
    float uy = 0.0;
    float uz = 0.0;
    float Ex = 0.0;
    float Ey = 0.0;
    float Ez = 0.0;
    float Bx = 0.0;
    float By = 0.0;
    float Bz = 0.0;

    fluidCell_3D_ideal_with_EM *HydroCell_3D_ideal_ptr1;
    fluidCell_3D_ideal_with_EM *HydroCell_3D_ideal_ptr2;
    for (int iptau = 0; iptau < 2; iptau++) {
        float taufactor;
        if (iptau == 0)
            taufactor = 1. - taufrac;
        else
            taufactor = taufrac;
        for (int ipeta = 0; ipeta < 2; ipeta++) {
            float etafactor;
            if (ipeta == 0)
                etafactor = 1. - etafrac;
            else
                etafactor = etafrac;
            for (int ipy = 0; ipy < 2; ipy++) {
                float yfactor;
                if (ipy == 0)
                    yfactor = 1. - yfrac;
                else
                    yfactor = yfrac;

                float prefrac = yfactor*etafactor*taufactor;

                HydroCell_3D_ideal_ptr1 = (
                    &lattice_3D_ideal_withEM_[idx_map_[position[0][ipy][ipeta][iptau]]]);
                HydroCell_3D_ideal_ptr2 = (
                    &lattice_3D_ideal_withEM_[idx_map_[position[1][ipy][ipeta][iptau]]]);
                T += prefrac*(
                    (1. - xfrac)*HydroCell_3D_ideal_ptr1->temperature
                    + xfrac*HydroCell_3D_ideal_ptr2->temperature);
                ed += prefrac*(
                    (1. - xfrac)*HydroCell_3D_ideal_ptr1->ed
                    + xfrac*HydroCell_3D_ideal_ptr2->ed);
                p += prefrac*(
                    (1. - xfrac)*HydroCell_3D_ideal_ptr1->pressure
                    + xfrac*HydroCell_3D_ideal_ptr2->pressure);
                ux += prefrac*((1. - xfrac)*HydroCell_3D_ideal_ptr1->ux
                                + xfrac*HydroCell_3D_ideal_ptr2->ux);
                uy += prefrac*((1. - xfrac)*HydroCell_3D_ideal_ptr1->uy
                                + xfrac*HydroCell_3D_ideal_ptr2->uy);
                uz += prefrac*((1. - xfrac)*HydroCell_3D_ideal_ptr1->uz
                               + xfrac*HydroCell_3D_ideal_ptr2->uz);
                Ex += prefrac*((1. - xfrac)*HydroCell_3D_ideal_ptr1->Ex
                                + xfrac*HydroCell_3D_ideal_ptr2->Ex);
                Ey += prefrac*((1. - xfrac)*HydroCell_3D_ideal_ptr1->Ey
                                + xfrac*HydroCell_3D_ideal_ptr2->Ey);
                Ez += prefrac*((1. - xfrac)*HydroCell_3D_ideal_ptr1->Ez
                               + xfrac*HydroCell_3D_ideal_ptr2->Ez);
                Bx += prefrac*((1. - xfrac)*HydroCell_3D_ideal_ptr1->Bx
                                + xfrac*HydroCell_3D_ideal_ptr2->Bx);
                By += prefrac*((1. - xfrac)*HydroCell_3D_ideal_ptr1->By
                                + xfrac*HydroCell_3D_ideal_ptr2->By);
                Bz += prefrac*((1. - xfrac)*HydroCell_3D_ideal_ptr1->Bz
                               + xfrac*HydroCell_3D_ideal_ptr2->Bz);
            }
        }
    }

    float ut = sqrt(1. + ux*ux + uy*uy + uz*uz);
    vx = ux/ut;
    vy = uy/ut;
    vz = uz/ut;

    info->temperature = T;
    info->vx = vx;
    info->vy = vy;
    info->vz = vz;

    info->ed = ed;
    info->sd = (ed + p)/(T + 1e-16);
    info->pressure = p;

    info->Ex = Ex;
    info->Ey = Ey;
    info->Ez = Ez;

    info->Bx = Bx;
    info->By = By;
    info->Bz = Bz;

    return;
}


void Hydroinfo_MUSIC::output_temperature_evolution(string filename_base) {
    fluidCell *hydroInfo = new fluidCell;
    for (int i = 0; i < itaumax; i++) {
        float tau = hydroTau0 + i*hydroDtau;
        ostringstream filename;
        filename << filename_base << "_tau_" << tau << ".dat";
        ofstream temp_evo(filename.str().c_str());
        for (int ix = 0; ix < ixmax; ix++) {
            float x_local = -hydroXmax + ix*hydroDx;
            for (int iy = 0; iy < ixmax; iy++) {
                float y_local = -hydroXmax + iy*hydroDx;
                getHydroValues(x_local, y_local, 0.0, tau, hydroInfo);
                float temp_local = hydroInfo->temperature;
                temp_evo << scientific << setw(16) << setprecision(8)
                         << temp_local << "   ";
            }
            temp_evo << endl;
        }
        temp_evo.close();
    }
    delete hydroInfo;
}


void Hydroinfo_MUSIC::update_grid_info(
    float tau0, float tau_max, float dtau,
    float x_max, float dx, float eta_max, float deta) {
    hydroTau0 = tau0;
    hydroTauMax = tau_max;
    hydroDtau = dtau;
    hydroXmax = x_max;
    hydroDx = dx;
    hydro_eta_max = eta_max;
    hydroDeta = deta;
    if (hydroWhichHydro == 8) {
        itaumax = static_cast<int>((tau_max-tau0)/dtau+0.001);
        ixmax = static_cast<int>(2*x_max/dx+0.001);
        ietamax = static_cast<int>(2*eta_max/deta+0.001);
    }
    if (hydroWhichHydro == 6) {
        itaumax = static_cast<int>((tau_max-tau0)/dtau+0.001);
        ixmax = static_cast<int>(2*x_max/dx+0.001);
        ietamax = static_cast<int>(2*eta_max/deta+0.001);
    }
}
