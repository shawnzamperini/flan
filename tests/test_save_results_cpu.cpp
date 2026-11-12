/*
* Tests for saving results to netCDF (save_results.cpp): a single-zone
* background keeps the original layout, and a multi-zone (SOLEDGE3X)
* background gets a zones group and one group per zone, with interior
* cells only. A sputtered source is saved in a sputtering group.
*/

#include <cmath>
#include <filesystem>
#include <gtest/gtest.h>
#include <string>
#include <vector>

#include "background.h"
#include "boundary.h"
#include "impurity_stats.h"
#include "netcdf"
#include "options.h"
#include "read_bkg.h"
#include "read_s3xe.h"
#include "s3x_test_files.h"
#include "save_results.h"
#include "sputtering.h"


namespace
{
	class SaveResultsTest : public testing::Test
	{
	protected:
		std::filesystem::path m_dir {std::filesystem::temp_directory_path()
			/ ("flan_test_save_" + std::string(testing::UnitTest::
			GetInstance()->current_test_info()->name()))};

		void SetUp() override {std::filesystem::create_directories(m_dir);}
		void TearDown() override {std::filesystem::remove_all(m_dir);}

		// Path of the saved file for a case name
		std::string nc_path() const {return (m_dir / "case.nc").string();}
		std::string case_name() const {return (m_dir / "case").string();}
	};

	// Reads a whole variable
	template <typename T>
	std::vector<T> read_var(const netCDF::NcGroup& group,
		const std::string& name)
	{
		const netCDF::NcVar var {group.getVar(name)};
		std::size_t n {1};
		for (const auto& d : var.getDims()) n *= d.getSize();
		std::vector<T> out(n);
		var.getVar(out.data());
		return out;
	}
}

// A single-zone background keeps its grid dimensions and groups at the
// top level, with no zones group
TEST_F(SaveResultsTest, SingleZoneLayout)
{
	Options::Options opts {};
	opts.set_bkg_source("test");
	opts.set_test_opt("gradb");
	opts.set_case_name(case_name());
	const auto bkg {Background::read_bkg(opts)};
	Impurity::Statistics stats {bkg};
	stats.calc_density(bkg, 1);
	SaveResults::save_results(bkg, stats, opts);

	const netCDF::NcFile file {nc_path(), netCDF::NcFile::read};
	EXPECT_FALSE(file.getDim("x").isNull());
	EXPECT_FALSE(file.getDim("grid z").isNull());
	EXPECT_EQ(file.getDim("x").getSize(),
		static_cast<std::size_t>(bkg.get_dim2()));
	EXPECT_TRUE(file.getGroup("zones").isNull());
	EXPECT_FALSE(file.getGroup("geometry").getVar("x").isNull());
	EXPECT_FALSE(file.getGroup("output").getVar("Nz").isNull());
}

// A multi-zone background gets a zones group with the face table, and a
// group per zone holding that zone's interior values
TEST_F(SaveResultsTest, MultiZoneLayout)
{
	// Write and read the small SOLEDGE3X case
	const S3XTest::TestMesh mesh {};
	S3XTest::write_case(m_dir.string(), mesh,
		[&](const std::string& path)
		{return S3X::read_mesh(path, mesh.nphi).zones;});
	Options::Options opts {};
	opts.set_bkg_source("soledge3x");
	opts.set_s3x_dir(m_dir.string());
	opts.set_s3x_frame_start(0);
	opts.set_s3x_frame_end(1);
	opts.set_case_name(case_name());
	const auto bkg {Background::read_bkg(opts)};

	// Statistics with a recognizable count and weight in zone 2
	Impurity::Statistics stats {bkg};
	stats.get_counts()(1, 1, 4, 7, 2) = 5;
	stats.get_weights()(1, 1, 4, 7, 2) = 3.0;
	stats.calc_density(bkg, 1);
	SaveResults::save_results(bkg, stats, opts);

	const netCDF::NcFile file {nc_path(), netCDF::NcFile::read};
	EXPECT_TRUE(file.getDim("x").isNull());
	EXPECT_EQ(file.getDim("time").getSize(), 2u);
	EXPECT_EQ(read_var<double>(file.getGroup("geometry"), "time"),
		bkg.get_times());
	EXPECT_FALSE(file.getGroup("input").getVar("s3x_dir").isNull());

	// Zone table: sizes, and faces numbered from 1 with -1 absorbing, -2
	// core and -4 wall, rows xmin, xmax, ymin, ymax, zmin, zmax
	const netCDF::NcGroup zones_group {file.getGroup("zones")};
	ASSERT_FALSE(zones_group.isNull());
	EXPECT_EQ(read_var<int>(zones_group, "nx"),
		(std::vector<int> {12, 12, 8}));
	EXPECT_EQ(read_var<int>(zones_group, "ny"),
		(std::vector<int> {40, 40, 40}));
	EXPECT_EQ(read_var<int>(zones_group, "faces"), (std::vector<int> {
		-2, -2, 1,
		3, -1, -1,
		2, 1, -4,
		2, 1, -4,
		1, 2, 3,
		1, 2, 3}));

	// Each zone has its own dimensions and groups
	for (int zone {}; zone < 3; ++zone)
	{
		const netCDF::NcGroup g {file.getGroup("zone"
			+ std::to_string(zone + 1))};
		ASSERT_FALSE(g.isNull()) << zone;
		EXPECT_EQ(g.getDim("x").getSize(),
			static_cast<std::size_t>(bkg.get_zones().get_nx(zone)));
		EXPECT_EQ(g.getDim("grid y").getSize(), 41u);
		for (const char* sub : {"geometry", "background", "output"})
			EXPECT_FALSE(g.getGroup(sub).isNull()) << zone << " " << sub;
	}

	// Values are the zone's interior cells, in (t, x, y, z) order
	const netCDF::NcGroup z2 {file.getGroup("zone2")};
	const int nx {12}, ny {40}, nz {mesh.nphi};
	auto idx4 = [&](int t, int i, int j, int k)
		{return ((t * nx + i) * ny + j) * nz + k;};
	auto idx3 = [&](int i, int j, int k) {return (i * ny + j) * nz + k;};

	const auto ne {read_var<double>(z2.getGroup("background"), "ne")};
	ASSERT_EQ(ne.size(), 2u * nx * ny * nz);
	EXPECT_DOUBLE_EQ(ne[idx4(0, 0, 0, 0)], bkg.get_ne()(1, 0, 0, 0, 0));
	EXPECT_DOUBLE_EQ(ne[idx4(1, 11, 39, 5)], bkg.get_ne()(1, 1, 11, 39, 5));

	const auto X {read_var<double>(z2.getGroup("geometry"), "X")};
	EXPECT_DOUBLE_EQ(X[idx3(3, 17, 2)], bkg.get_X()(1, 3, 17, 2));
	const auto chi {read_var<double>(z2.getGroup("geometry"), "chi")};
	EXPECT_EQ(chi.size(), static_cast<std::size_t>(nx * ny * nz));

	const auto Nz {read_var<int>(z2.getGroup("output"), "Nz")};
	EXPECT_EQ(Nz[idx4(1, 4, 7, 2)], 5);
	const auto nz_out {read_var<double>(z2.getGroup("output"), "nz")};
	EXPECT_DOUBLE_EQ(nz_out[idx4(1, 4, 7, 2)],
		stats.get_density()(1, 1, 4, 7, 2));

	// Data the background doesn't have is left out
	EXPECT_TRUE(z2.getGroup("geometry").getVar("gij_00").isNull());
	EXPECT_TRUE(z2.getGroup("background").getVar("B_R").isNull());

	// Poloidal corners: the first corner of zone 3 is at r = r_mid,
	// theta = 0, so R = R0 + a0 r_mid and Z = 0
	const netCDF::NcGroup geo3 {file.getGroup("zone3").getGroup("geometry")};
	const auto R {read_var<double>(geo3, "R_corners")};
	const auto Z {read_var<double>(geo3, "Z_corners")};
	ASSERT_EQ(R.size(), 9u * 41u);
	EXPECT_NEAR(R[0], mesh.R0 + mesh.a0 * mesh.r_mid, 1e-12);
	EXPECT_NEAR(Z[0], 0.0, 1e-12);
}

// The sputtered source of a SOLEDGE3X background is saved in a sputtering
// group, with an element dimension, per-element values and per frame and
// element values
TEST_F(SaveResultsTest, Sputtering)
{
	const S3XTest::TestMesh mesh {};
	S3XTest::write_case(m_dir.string(), mesh,
		[&](const std::string& path)
		{return S3X::read_mesh(path, mesh.nphi).zones;});
	Options::Options opts {};
	opts.set_bkg_source("soledge3x");
	opts.set_s3x_dir(m_dir.string());
	opts.set_s3x_frame_start(0);
	opts.set_s3x_frame_end(1);
	opts.set_case_name(case_name());
	opts.set_imp_source_opt("sputter");
	opts.set_sput_table_file(std::string {FLAN_SOURCE_DIR}
		+ "/data/sputtering/D_on_W.h5");
	const auto bkg {Background::read_bkg(opts)};

	const Sputtering::Surface surface(bkg);
	const Sputtering::YieldTable table(opts.sput_table_file());
	const Sputtering::Source source(surface, bkg, table, opts);
	Sputtering::Deposition dep(source.surface(), 74);
	dep.add_into_wall_cell(source.surface().elements()[3].wall_zone,
		source.surface().elements()[3].wall_i,
		source.surface().elements()[3].wall_j,
		source.surface().elements()[3].k, 2, source.surface().elements()[3].i,
		source.surface().elements()[3].j, 1e-5, 2, 7.0);
	Impurity::Statistics stats {bkg};
	stats.calc_density(bkg, 1);
	Boundary::Fates fates {};
	fates.chi = 151;
	fates.max_time = 47;
	SaveResults::save_results(bkg, stats, opts, &source, &dep, &fates);

	const netCDF::NcFile file {nc_path(), netCDF::NcFile::read};
	const netCDF::NcGroup group {file.getGroup("sputtering")};
	ASSERT_FALSE(group.isNull());
	const std::size_t n {source.size()};
	ASSERT_EQ(group.getDim("element").getSize(), n);

	// Zone 3's last psi row is the wall, seen from the row below it
	for (const int zone : read_var<int>(group, "zone")) EXPECT_EQ(zone, 3);
	EXPECT_EQ(read_var<int>(group, "i")[0], mesh.npsi_sol - 2);
	EXPECT_EQ(read_var<double>(group, "area")[5],
		surface.elements()[5].area);

	// Per frame and element, frame by frame
	const auto rate {read_var<double>(group, "rate")};
	ASSERT_EQ(rate.size(), 2 * n);
	EXPECT_EQ(rate, source.rate());
	EXPECT_EQ(read_var<double>(group, "flux_exb"), source.flux_exb());
	EXPECT_EQ(read_var<double>(group, "total_rate"), source.total_rate());
	EXPECT_GT(source.mean_total_rate(), 0.0);

	std::string table_file {};
	group.getAtt("sput_table_file").getValues(table_file);
	EXPECT_EQ(table_file, opts.sput_table_file());

	// Where the atoms landed: the source averaged over the frames, what
	// landed on each element, and the difference
	const auto mean {read_var<double>(group, "source_mean")};
	const auto redep {read_var<double>(group, "redeposited")};
	const auto net {read_var<double>(group, "net_erosion")};
	EXPECT_EQ(mean, source.mean_rate());
	EXPECT_DOUBLE_EQ(redep[3], 7.0);
	EXPECT_DOUBLE_EQ(net[3], mean[3] - 7.0);
	EXPECT_EQ(read_var<double>(group, "redeposited_charge")[2], 7.0);
	double total {};
	group.getAtt("redeposited_total").getValues(&total);
	EXPECT_DOUBLE_EQ(total, 7.0);

	// What happened to the particles, as counts
	long long chi {}, max_time {}, core {-1};
	group.getAtt("fates_chi").getValues(&chi);
	group.getAtt("fates_max_time").getValues(&max_time);
	group.getAtt("fates_core").getValues(&core);
	EXPECT_EQ(chi, 151);
	EXPECT_EQ(max_time, 47);
	EXPECT_EQ(core, 0);

	// No sputtering group without a source
	SaveResults::save_results(bkg, stats, opts);
	const netCDF::NcFile again {nc_path(), netCDF::NcFile::read};
	EXPECT_TRUE(again.getGroup("sputtering").isNull());
}
