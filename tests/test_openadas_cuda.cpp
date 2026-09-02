#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>

#include "background.h"
#include "background_device.h"
#include "options.h"
#include "openadas.h"
#include "openadas_device.h"
#include "pcg32.h"
#include "read_bkg.h"
#include "slots.h"
#include "slots_device.h"

#include "init_rngs.cuh"
#include "openadas.cuh"


namespace fs = std::filesystem;

// Write a minimal valid multi-charge ADAS file for Helium (Z=2)
static void write_fake_adas_helium(const std::string& path,
                                   double rate00_c0,
                                   double rate01_c0,
                                   double rate10_c0,
                                   double rate11_c0,
                                   double rate00_c1,
                                   double rate01_c1,
                                   double rate10_c1,
                                   double rate11_c1)
{
    std::ofstream f(path);

    // atomic_number=2, ndens=2, ntemp=2, charge_low=0, charge_high=1
    f << "  2   2   2    0    1     /HELIUM             / TEST\n";
    f << " -------------------------------------------------------------------------------\n";

    // Te grid (2 values)
    f << "  10.0  20.0\n";

    // Ne grid (2 values)
    f << "   0.1   1.0\n";

    // Charge state 0 block (He0 → He1+ for SCD, He1+ → He0 for ACD)
    f << " ---------------------------------------------------/ Z0= 0   / DATE= 05.04.90\n";
    f << " " << rate00_c0 << " " << rate01_c0 << "\n";
    f << " " << rate10_c0 << " " << rate11_c0 << "\n";

    // Charge state 1 block (He1+ → He2+ for SCD, He2+ → He1+ for ACD)
    f << " ---------------------------------------------------/ Z1= 1   / DATE= 05.04.90\n";
    f << " " << rate00_c1 << " " << rate01_c1 << "\n";
    f << " " << rate10_c1 << " " << rate11_c1 << "\n";

    f.close();
}

static void setup_single_particle(Slots::Slots& slots,
    const Background::Background& bkg, int init_q)
{
    slots.set_q(0, init_q);
    slots.set_Z(2);  // Helium (Z=2)

    slots.set_t(0, bkg.get_times()[0]);
    slots.set_x(0, 0.01);
    slots.set_y(0, 0.0);
    slots.set_z(0, 0.0);

    slots.set_tidx(0, 0);
    slots.set_xidx(0, 10);
    slots.set_yidx(0, 10);
    slots.set_zidx(0, 4);
}

TEST(IonizRecomb, IonizationOnly_IncreasesCharge_He)
{
    fs::path temp_dir = fs::temp_directory_path() / "fake_adas_ioniz_he";
    fs::create_directories(temp_dir / "scd2");
    fs::create_directories(temp_dir / "acd2");

    // Ionization: large rates for both charge states
    write_fake_adas_helium((temp_dir / "scd2/scd2_he.dat").string(),
                           -1.0, -1.0,
                           -1.0, -1.0,
                           -1.0, -1.0,
                           -1.0, -1.0);

    // Recombination: tiny rates for both charge states
    write_fake_adas_helium((temp_dir / "acd2/acd2_he.dat").string(),
                           -99.0, -99.0,
                           -99.0, -99.0,
                           -99.0, -99.0,
                           -99.0, -99.0);

    OpenADAS::OpenADAS oa_ioniz(temp_dir.string(), 2, 2, "scd");
    OpenADAS::OpenADAS oa_recomb(temp_dir.string(), 2, 2, "acd");

	// Copy to device
	OpenADAS::OpenADASDevice oa_ioniz_d {oa_ioniz.to_device(0)};
	OpenADAS::OpenADASDevice oa_recomb_d {oa_recomb.to_device(0)};

    Options::Options opts {};
    opts.set_bkg_source("test");
    opts.set_test_opt("gyrate");

	// Create background and copy to device
    Background::Background bkg = Background::read_bkg(opts);
	Background::BackgroundDevice bkg_d {bkg.to_device()};

    // Start with neutral helium (He0)
    Slots::Slots slots(1);
    setup_single_particle(slots, bkg, 0);

	// Copy to device
	Slots::SlotsDevice slots_d {slots.to_device()};

	// Setup the RNGs
	pcg32* rngs_d {init_rngs_cuda(slots_d, 4)};

    int ioniz_warnings = 0;
    int recomb_warnings = 0;

	OpenADAS::ioniz_recomb_gpu(slots_d, bkg_d, oa_ioniz_d, 
		oa_recomb_d, 1e-6, ioniz_warnings, recomb_warnings, rngs_d);
	
	// Copy back to host, copies slots_d into slots. Free memory on GPU.
	slots = slots.to_host(slots_d);
	Slots::free_slots(slots_d);
	OpenADAS::free_oa(oa_ioniz_d, 0);
	OpenADAS::free_oa(oa_recomb_d, 0);

    EXPECT_EQ(slots.q()[0], 1);   // He0 → He1+

	// Should be 1 ionization warning since we've made the ionization rates
	// artifically large
    EXPECT_EQ(ioniz_warnings, 1);
    EXPECT_EQ(recomb_warnings, 0);
}


TEST(IonizRecomb, RecombinationOnly_DecreasesCharge_He)
{
    fs::path temp_dir = fs::temp_directory_path() / "fake_adas_recomb_he";
    fs::create_directories(temp_dir / "scd2");
    fs::create_directories(temp_dir / "acd2");

    // Ionization: tiny rates
    write_fake_adas_helium((temp_dir / "scd2/scd2_he.dat").string(),
                           -99.0, -99.0,
                           -99.0, -99.0,
                           -99.0, -99.0,
                           -99.0, -99.0);

    // Recombination: large rates
    write_fake_adas_helium((temp_dir / "acd2/acd2_he.dat").string(),
                           -1.0, -1.0,
                           -1.0, -1.0,
                           -1.0, -1.0,
                           -1.0, -1.0);

    OpenADAS::OpenADAS oa_ioniz(temp_dir.string(), 2, 2, "scd");
    OpenADAS::OpenADAS oa_recomb(temp_dir.string(), 2, 2, "acd");

	// Copy to device
	OpenADAS::OpenADASDevice oa_ioniz_d {oa_ioniz.to_device(0)};
	OpenADAS::OpenADASDevice oa_recomb_d {oa_recomb.to_device(0)};

    Options::Options opts {};
    opts.set_bkg_source("test");
    opts.set_test_opt("gyrate");

	// Create background and copy to device
    Background::Background bkg = Background::read_bkg(opts);
	Background::BackgroundDevice bkg_d {bkg.to_device()};

    // Start with He+
    Slots::Slots slots(1);
    setup_single_particle(slots, bkg, 1);

	// Copy to device
	Slots::SlotsDevice slots_d {slots.to_device()};

	// Setup the RNGs
	pcg32* rngs_d {init_rngs_cuda(slots_d, 4)};

    int ioniz_warnings = 0;
    int recomb_warnings = 0;

	OpenADAS::ioniz_recomb_gpu(slots_d, bkg_d, oa_ioniz_d, 
		oa_recomb_d, 1e-6, ioniz_warnings, recomb_warnings, rngs_d);
	
	// Copy back to host, copies slots_d into slots. Free memory on GPU.
	slots = slots.to_host(slots_d);
	Slots::free_slots(slots_d);
	OpenADAS::free_oa(oa_ioniz_d, 0);
	OpenADAS::free_oa(oa_recomb_d, 0);

    EXPECT_EQ(slots.q()[0], 0);   // He1+ → He0
    EXPECT_EQ(ioniz_warnings, 0);

	// Similar to previous test, should be a recombination warning since we
	// made the recombination rates arbitrarily large
    EXPECT_EQ(recomb_warnings, 1);
}

TEST(IonizRecomb, CancellationCase_NoNetChange_He)
{
    fs::path temp_dir = fs::temp_directory_path() / "fake_adas_cancel_he";
    fs::create_directories(temp_dir / "scd2");
    fs::create_directories(temp_dir / "acd2");

    // Ionization: large rates
    write_fake_adas_helium((temp_dir / "scd2/scd2_he.dat").string(),
                           -1.0, -1.0,
                           -1.0, -1.0,
                           -1.0, -1.0,
                           -1.0, -1.0);

    // Recombination: large rates
    write_fake_adas_helium((temp_dir / "acd2/acd2_he.dat").string(),
                           -1.0, -1.0,
                           -1.0, -1.0,
                           -1.0, -1.0,
                           -1.0, -1.0);

    OpenADAS::OpenADAS oa_ioniz(temp_dir.string(), 2, 2, "scd");
    OpenADAS::OpenADAS oa_recomb(temp_dir.string(), 2, 2, "acd");

	// Copy to device
	OpenADAS::OpenADASDevice oa_ioniz_d {oa_ioniz.to_device(0)};
	OpenADAS::OpenADASDevice oa_recomb_d {oa_recomb.to_device(0)};

    Options::Options opts {};
    opts.set_bkg_source("test");
    opts.set_test_opt("gyrate");

	// Create background and copy to device
    Background::Background bkg = Background::read_bkg(opts);
	Background::BackgroundDevice bkg_d {bkg.to_device()};

    // Start with He1+
    Slots::Slots slots(1);
    setup_single_particle(slots, bkg, 1);

	// Copy to device
	Slots::SlotsDevice slots_d {slots.to_device()};

	// Setup the RNGs
	pcg32* rngs_d {init_rngs_cuda(slots_d, 4)};

    int ioniz_warnings = 0;
    int recomb_warnings = 0;

	OpenADAS::ioniz_recomb_gpu(slots_d, bkg_d, oa_ioniz_d, 
		oa_recomb_d, 1e-6, ioniz_warnings, recomb_warnings, rngs_d);

	// Copy back to host, copies slots_d into slots. Free memory on GPU.
	slots = slots.to_host(slots_d);
	Slots::free_slots(slots_d);
	OpenADAS::free_oa(oa_ioniz_d, 0);
	OpenADAS::free_oa(oa_recomb_d, 0);

    EXPECT_EQ(slots.q()[0], 1);   // dq = +1 - 1 = 0
    EXPECT_EQ(ioniz_warnings, 1);
    EXPECT_EQ(recomb_warnings, 1);
}
