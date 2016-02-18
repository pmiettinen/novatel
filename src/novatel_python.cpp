/*
* Copyright (C) 2016 Swift Navigation Inc.
* Contact: Pasi Miettinen <pasi.miettinen@exafore.com>
*
* This source is subject to the license found in the file 'LICENSE'
* which must be be distributed together with this source. All other
* rights reserved.
*
* THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
* KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR
* PURPOSE.
*/

#include "novatel/novatel_python.h"

using namespace novatel;

boost::posix_time::time_duration::tick_type MillisecondsSinceEpoch()
{
  return (boost::posix_time::microsec_clock::universal_time() -\
    boost::posix_time::ptime(boost::gregorian::date(1980, 1, 6))).\
    total_milliseconds();
}

NovatelPython::NovatelPython(std::string id)
  : Novatel()
  , m_id(id)
{
  set_time_handler(MillisecondsSinceEpoch);
}

NovatelPython::NovatelPython(std::string id, std::string raw_log_name)
  : Novatel()
  , m_id(id)
{
  set_time_handler(MillisecondsSinceEpoch);
  CreateRawLog(raw_log_name);
}


NovatelPython::~NovatelPython()
{
}

PyObject *
NovatelPython::python_set_best_position_callback(PyObject *callback)
{
  python_best_position_callback_ = callback;
  set_best_position_callback(
    boost::bind(&NovatelPython::best_position_callback, this, _1, _2));
  return Py_None;
}

PyObject *
NovatelPython::python_set_best_position_ecef_callback(
  PyObject *callback)
{
  python_best_position_ecef_callback_ = callback;
  set_best_position_ecef_callback(
    boost::bind(&NovatelPython::best_position_ecef_callback, this, _1,
    _2));
  return Py_None;
}

PyObject *
NovatelPython::python_set_range_measurements_callback(
  PyObject *callback)
{
  python_range_measurements_callback_ = callback;
  set_range_measurements_callback(
    boost::bind(&NovatelPython::range_measurements_callback, this, _1,
    _2));
  return Py_None;
}

PyObject *
NovatelPython::python_set_raw_gps_word_callback(PyObject *callback)
{
  python_raw_gps_word_callback_ = callback;
  set_raw_gps_word_callback(
    boost::bind(&NovatelPython::raw_gps_word_callback, this, _1, _2));
  return Py_None;
}

PyObject *
NovatelPython::python_set_gps_ephem_callback(PyObject *callback)
{
  python_gps_ephem_callback_ = callback;
  set_gps_ephemeris_callback(
    boost::bind(&NovatelPython::gps_ephem_callback, this, _1, _2));
  return Py_None;
}

PyObject *
NovatelPython::python_set_gps_almanac_callback(PyObject *callback)
{
  python_gps_almanac_callback_ = callback;
  set_almanac_callback(
    boost::bind(&NovatelPython::gps_almanac_callback, this, _1, _2));
  return Py_None;
}

PyObject *
NovatelPython::python_set_best_sats_callback(PyObject *callback)
{
  python_best_sats_callback_ = callback;
  set_best_sats_callback(
    boost::bind(&NovatelPython::best_sats_callback, this, _1, _2));
  return Py_None;
}

void NovatelPython::best_position_callback(Position const &best_pos,
  double const &timestamp)
{
  call_python(python_best_position_callback_, PositionPython(best_pos),
    timestamp);
}

void NovatelPython::best_position_ecef_callback(
  PositionEcef const &best_pos_ecef, double const &timestamp)
{
  call_python(python_best_position_ecef_callback_,
    PositionEcefPython(best_pos_ecef), timestamp);
}

void NovatelPython::range_measurements_callback(
  RangeMeasurements const &range_meas, double const &timestamp)
{
  call_python(python_range_measurements_callback_,
    RangeMeasurementsPython(range_meas), timestamp);
}

void NovatelPython::raw_gps_word_callback(
  RawGpsWord const &raw_gps_word, double const &timestamp)
{
  call_python(python_raw_gps_word_callback_,
    RawGpsWordPython(raw_gps_word), timestamp);
}

void NovatelPython::gps_ephem_callback(
  GpsEphemeris const &gps_ephem, double const &timestamp)
{
  call_python(python_gps_ephem_callback_, GpsEphemerisPython(gps_ephem),
    timestamp);
}

void NovatelPython::gps_almanac_callback(
  Almanac const &gps_almanac, double const &timestamp)
{
  call_python(python_gps_almanac_callback_, AlmanacPython(gps_almanac),
    timestamp);
}

void NovatelPython::best_sats_callback(
  BestSats const &best_sats, double const &timestamp)
{
  call_python(python_best_sats_callback_, BestSatsPython(best_sats),
    timestamp);
}

template<typename T>
void NovatelPython::call_python(PyObject *callable, T const &data,
  double const &timestamp)
{
  if (callable) {
    // Ensure that the current thread is ready to call the Python C API
    PyGILState_STATE state = PyGILState_Ensure();

    // invoke the python function
    try {
      boost::python::call<void>(callable, boost::ref(data), timestamp,
        m_id);
    }
    catch(error_already_set &) {
      PyObject *ptype, *pvalue, *ptraceback;
      PyErr_Fetch(&ptype, &pvalue, &ptraceback);
      //pvalue contains error message
      //ptraceback contains stack snapshot and many other information
      //(see python traceback structure)

      //Get error message
      std::cout << PyString_AsString(pvalue);
    }

    // release the global interpreter lock
    // so other threads can resume execution
    PyGILState_Release(state);
  }
}

BOOST_PYTHON_MODULE(novatel)
{
  class_<Novatel, boost::noncopyable>("Novatel")
    .def("IsConnected", &Novatel::IsConnected)
    .def("Connect", &Novatel::Connect)
    .def("Disconnect", &Novatel::Disconnect)
    .def("UnlogAll", &Novatel::UnlogAll)
    .def("ConfigureLogs", &Novatel::ConfigureLogs)
  ;

  class_<NovatelPython, bases<Novatel>,
    boost::noncopyable>("NovatelPython", init<std::string>())
    .def(init<std::string, std::string>())
    .def("python_set_best_position_callback",
      &NovatelPython::python_set_best_position_callback)
    .def("python_set_best_position_ecef_callback",
      &NovatelPython::python_set_best_position_ecef_callback)
    .def("python_set_range_measurements_callback",
      &NovatelPython::python_set_range_measurements_callback)
    .def("python_set_raw_gps_word_callback",
      &NovatelPython::python_set_raw_gps_word_callback)
    .def("python_set_gps_ephem_callback",
      &NovatelPython::python_set_gps_ephem_callback)
    .def("python_set_gps_almanac_callback",
      &NovatelPython::python_set_gps_almanac_callback)
    .def("python_set_best_sats_callback",
      &NovatelPython::python_set_best_sats_callback)
  ;

  class_<Oem4BinaryHeader>("Oem4BinaryHeader")
    .def_readonly("sync1", &Oem4BinaryHeader::sync1)
    .def_readonly("sync2", &Oem4BinaryHeader::sync2)
    .def_readonly("sync3", &Oem4BinaryHeader::sync3)
  ;

  class_<Oem4BinaryHeaderPython>("Oem4BinaryHeaderPython",
    init<Oem4BinaryHeader>())
    .add_property("sync1", &Oem4BinaryHeaderPython::sync1)
    .add_property("sync2", &Oem4BinaryHeaderPython::sync2)
    .add_property("sync3", &Oem4BinaryHeaderPython::sync3)
    .add_property("message_id", &Oem4BinaryHeaderPython::message_id)
    .add_property("message_type", &Oem4BinaryHeaderPython::message_type)
    .add_property("port_address", &Oem4BinaryHeaderPython::port_address)
    .add_property("message_length",
      &Oem4BinaryHeaderPython::message_length)
    .add_property("sequence", &Oem4BinaryHeaderPython::sequence)
    .add_property("idle", &Oem4BinaryHeaderPython::idle)
    .add_property("time_status", &Oem4BinaryHeaderPython::time_status)
    .add_property("gps_week", &Oem4BinaryHeaderPython::gps_week)
    .add_property("gps_millisecs",
      &Oem4BinaryHeaderPython::gps_millisecs)
    .add_property("status", &Oem4BinaryHeaderPython::status)
    .add_property("Reserved", &Oem4BinaryHeaderPython::Reserved)
    .add_property("version", &Oem4BinaryHeaderPython::version)
  ;

  enum_<SolutionStatus>("SolutionStatus")
    .value("SOL_COMPUTED", SOL_COMPUTED)
    .value("INSUFFICIENT_OBS", INSUFFICIENT_OBS)
    .value("NO_CONVERGENCE", NO_CONVERGENCE)
    .value("SINGULARITY", SINGULARITY)
    .value("COV_TRACE", COV_TRACE)
    .value("TEST_DIST", TEST_DIST)
    .value("COLD_START", COLD_START)
    .value("V_H_LIMIT", V_H_LIMIT)
    .value("VARIANCE", VARIANCE)
    .value("RESIDUALS", RESIDUALS)
    .value("DELTA_POS", DELTA_POS)
    .value("NEGATIVE_VAR", NEGATIVE_VAR)
    .value("INTEGRITY_WARNING", INTEGRITY_WARNING)
    .value("INS_INACTIVE", INS_INACTIVE)
    .value("INS_ALIGNING", INS_ALIGNING)
    .value("INS_BAD", INS_BAD)
    .value("IMU_UNPLUGGED", IMU_UNPLUGGED)
    .value("PENDING", PENDING)
    .value("INVALID_FIX", INVALID_FIX)
    .value("UNAUTHORIZED", UNAUTHORIZED)
  ;

  enum_<PositionType>("PositionType")
    .value("NONE", NONE)
    .value("FIXEDPOS", FIXEDPOS)
    .value("FIXEDHEIGHT", FIXEDHEIGHT)
    .value("Reserved", Reserved)
    .value("FLOATCONV", FLOATCONV)
    .value("WIDELANE", WIDELANE)
    .value("NARROWLANE", NARROWLANE)
    .value("DOPPLER_VELOCITY", DOPPLER_VELOCITY)
    .value("SINGLE", SINGLE)
    .value("PSRDIFF", PSRDIFF)
    .value("WAAS", WAAS)
    .value("PROPOGATED", PROPOGATED)
    .value("OMNISTAR", OMNISTAR)
    .value("L1_FLOAT", L1_FLOAT)
    .value("IONOFREE_FLOAT", IONOFREE_FLOAT)
    .value("NARROW_FLOAT", NARROW_FLOAT)
    .value("L1_INT", L1_INT)
    .value("WIDE_INT", WIDE_INT)
    .value("NARROW_INT", NARROW_INT)
    .value("RTK_DIRECT_INS", RTK_DIRECT_INS)
    .value("INS", INS)
    .value("INS_PSRSP", INS_PSRSP)
    .value("INS_PSRDIFF", INS_PSRDIFF)
    .value("INS_RTKFLOAT", INS_RTKFLOAT)
    .value("INS_RTKFIXED", INS_RTKFIXED)
    .value("OMNISTAR_HP", OMNISTAR_HP)
    .value("OMNISTAR_XP", OMNISTAR_XP)
    .value("CDGPS", CDGPS)
  ;

  enum_<DatumID>("DatumID")
    .value("ADIND", ADIND)
    .value("ARC50", ARC50)
    .value("ARC60", ARC60)
    .value("AGD66", AGD66)
    .value("AGD84", AGD84)
    .value("BUKIT", BUKIT)
    .value("ASTRO", ASTRO)
    .value("CHATM", CHATM)
    .value("CARTH", CARTH)
    .value("CAPE", CAPE)
    .value("DJAKA", DJAKA)
    .value("EGYPT", EGYPT)
    .value("ED50", ED50)
    .value("ED79", ED79)
    .value("GUNSG", GUNSG)
    .value("GEO49", GEO49)
    .value("GRB36", GRB36)
    .value("GUAM", GUAM)
    .value("HAWAII", HAWAII)
    .value("KAUAI", KAUAI)
    .value("MAUI", MAUI)
    .value("OAHU", OAHU)
    .value("HERAT", HERAT)
    .value("HJORS", HJORS)
    .value("HONGK", HONGK)
    .value("HUTZU", HUTZU)
    .value("INDIA", INDIA)
    .value("IRE65", IRE65)
    .value("KERTA", KERTA)
    .value("KANDA", KANDA)
    .value("LIBER", LIBER)
    .value("LUZON", LUZON)
    .value("MINDA", MINDA)
    .value("MERCH", MERCH)
    .value("NAHR", NAHR)
    .value("NAD83", NAD83)
    .value("CANADA", CANADA)
    .value("ALASKA", ALASKA)
    .value("NAD27", NAD27)
    .value("CARIBB", CARIBB)
    .value("MEXICO", MEXICO)
    .value("CAMER", CAMER)
    .value("MINNA", MINNA)
    .value("OMAN", OMAN)
    .value("PUERTO", PUERTO)
    .value("QORNO", QORNO)
    .value("ROME", ROME)
    .value("CHUA", CHUA)
    .value("SAM56", SAM56)
    .value("SAM69", SAM69)
    .value("CAMPO", CAMPO)
    .value("SACOR", SACOR)
    .value("YACAR", YACAR)
    .value("TANAN", TANAN)
    .value("TIMBA", TIMBA)
    .value("TOKYO", TOKYO)
    .value("TRIST", TRIST)
    .value("VITI", VITI)
    .value("WAK60", WAK60)
    .value("WGS72", WGS72)
    .value("WGS84", WGS84)
    .value("ZANDE", ZANDE)
    .value("USER", USER)
    .value("CSRS", CSRS)
    .value("ADIM", ADIM)
    .value("ARSM", ARSM)
    .value("ENW", ENW)
    .value("HTN", HTN)
    .value("INDB", INDB)
    .value("INDI", INDI)
    .value("IRL", IRL)
    .value("LUZA", LUZA)
    .value("LUZB", LUZB)
    .value("NAHC", NAHC)
    .value("NASP", NASP)
    .value("OGBM", OGBM)
    .value("OHAA", OHAA)
    .value("OHAB", OHAB)
    .value("OHAC", OHAC)
    .value("OHAD", OHAD)
    .value("OHIA", OHIA)
    .value("OHIB", OHIB)
    .value("OHIC", OHIC)
    .value("OHID", OHID)
    .value("TIL", TIL)
    .value("TOYM", TOYM)
  ;

  enum_<true_false>("true_false")
    .value("FALSE", FALSE)
    .value("TRUE", TRUE)
  ;

  enum_<yes_no>("yes_no")
    .value("no", no)
    .value("yes", yes)
  ;

  class_<PositionPython>("PositionPython", init<Position>())
    .add_property("header", &PositionPython::header)
    .add_property("solution_status", &PositionPython::solution_status)
    .add_property("position_type", &PositionPython::position_type)
    .add_property("latitude", &PositionPython::latitude)
    .add_property("longitude", &PositionPython::longitude)
    .add_property("height", &PositionPython::height)
    .add_property("undulation", &PositionPython::undulation)
    .add_property("datum_id", &PositionPython::datum_id)
    .add_property("latitude_standard_deviation",
       &PositionPython::latitude_standard_deviation)
    .add_property("longitude_standard_deviation",
      &PositionPython::longitude_standard_deviation)
    .add_property("height_standard_deviation",
      &PositionPython::height_standard_deviation)
    .add_property("base_station_id", &PositionPython::base_station_id)
    .add_property("differential_age", &PositionPython::differential_age)
    .add_property("solution_age", &PositionPython::solution_age)
    .add_property("number_of_satellites",
      &PositionPython::number_of_satellites)
    .add_property("number_of_satellites_in_solution",
      &PositionPython::number_of_satellites_in_solution)
    .add_property("num_gps_plus_glonass_l1",
      &PositionPython::num_gps_plus_glonass_l1)
    .add_property("num_gps_plus_glonass_l2",
      &PositionPython::num_gps_plus_glonass_l2)
    .add_property("reserved", &PositionPython::reserved)
    .add_property("extended_solution_status",
      &PositionPython::extended_solution_status)
    .add_property("reserved2", &PositionPython::reserved2)
    .add_property("signals_used_mask",
      &PositionPython::signals_used_mask)
    .add_property("crc", &PositionPython::crc)
  ;

  class_<PositionEcefPython>("PositionEcef", init<PositionEcef>())
    .add_property("header", &PositionEcefPython::header)
    .add_property("solution_status",
      &PositionEcefPython::solution_status)
    .add_property("position_type", &PositionEcefPython::position_type)
    .add_property("x_position", &PositionEcefPython::x_position)
    .add_property("y_position", &PositionEcefPython::y_position)
    .add_property("z_position", &PositionEcefPython::z_position)
    .add_property("x_standard_deviation",
      &PositionEcefPython::x_standard_deviation)
    .add_property("y_standard_deviation",
      &PositionEcefPython::y_standard_deviation)
    .add_property("z_standard_deviation",
       &PositionEcefPython::z_standard_deviation)
    .add_property("velocity_status",
      &PositionEcefPython::velocity_status)
    .add_property("velocity_type", &PositionEcefPython::velocity_type)
    .add_property("x_velocity", &PositionEcefPython::x_velocity)
    .add_property("y_velocity", &PositionEcefPython::y_velocity)
    .add_property("z_velocity", &PositionEcefPython::z_velocity)
    .add_property("x_velocity_standard_deviation",
      &PositionEcefPython::x_velocity_standard_deviation)
    .add_property("y_velocity_standard_deviation",
      &PositionEcefPython::y_velocity_standard_deviation)
    .add_property("z_velocity_standard_deviation",
      &PositionEcefPython::z_velocity_standard_deviation)
    .add_property("base_station_id",
      &PositionEcefPython::base_station_id)
    .add_property("velocity_latency",
      &PositionEcefPython::velocity_latency)
    .add_property("differential_age",
      &PositionEcefPython::differential_age)
    .add_property("solution_age", &PositionEcefPython::solution_age)
    .add_property("number_of_satellites",
      &PositionEcefPython::number_of_satellites)
    .add_property("number_of_satellites_in_solution",
      &PositionEcefPython::number_of_satellites_in_solution)
    .add_property("reserved", &PositionEcefPython::reserved)
    .add_property("extended_solution_status",
      &PositionEcefPython::extended_solution_status)
    .add_property("reserved2", &PositionEcefPython::reserved2)
    .add_property("signals_used_mask",
      &PositionEcefPython::signals_used_mask)
    .add_property("crc", &PositionEcefPython::crc)
  ;

  class_<ChannelStatusPython>("ChannelStatusPython",
    init<ChannelStatus>())
    .add_property("tracking_state",
      &ChannelStatusPython::tracking_state)
    .add_property("sv_chan_num", &ChannelStatusPython::sv_chan_num)
    .add_property("phase_lock_flag",
      &ChannelStatusPython::phase_lock_flag)
    .add_property("parity_known_flag",
      &ChannelStatusPython::parity_known_flag)
    .add_property("code_locked_flag",
      &ChannelStatusPython::code_locked_flag)
    .add_property("correlator_type",
      &ChannelStatusPython::correlator_type)
    .add_property("satellite_sys", &ChannelStatusPython::satellite_sys)
    .add_property("reserved1", &ChannelStatusPython::reserved1)
    .add_property("grouping", &ChannelStatusPython::grouping)
    .add_property("signal_type", &ChannelStatusPython::signal_type)
    .add_property("forward_err_correction",
      &ChannelStatusPython::forward_err_correction)
    .add_property("primary_L1_chan",
      &ChannelStatusPython::primary_L1_chan)
    .add_property("carrier_phase_meas",
      &ChannelStatusPython::carrier_phase_meas)
    .add_property("reserved2", &ChannelStatusPython::reserved2)
    .add_property("prn_lock_flag", &ChannelStatusPython::prn_lock_flag)
    .add_property("channel_assignment",
      &ChannelStatusPython::channel_assignment)
  ;

  class_<RangeDataPython>("RangeDataPython", init<RangeData>())
    .add_property("satellite_prn", &RangeDataPython::satellite_prn)
    .add_property("glonass_frequency",
      &RangeDataPython::glonass_frequency)
    .add_property("pseudorange", &RangeDataPython::pseudorange)
    .add_property("pseudorange_standard_deviation",
      &RangeDataPython::pseudorange_standard_deviation)
    .add_property("accumulated_doppler",
      &RangeDataPython::accumulated_doppler)
    .add_property("accumulated_doppler_std_deviation",
      &RangeDataPython::accumulated_doppler_std_deviation)
    .add_property("doppler", &RangeDataPython::doppler)
    .add_property("carrier_to_noise",
      &RangeDataPython::carrier_to_noise)
    .add_property("locktime", &RangeDataPython::locktime)
    .add_property("channel_status_packed",
      &RangeDataPython::channel_status_packed)
    .add_property("channel_status", &RangeDataPython::channel_status)
  ;

  class_<RangeMeasurementsPython>("RangeMeasurementsPython",
    init<RangeMeasurements>())
    .add_property("header", &RangeMeasurementsPython::header)
    .add_property("number_of_observations",
      &RangeMeasurementsPython::number_of_observations)
    .add_property("range_data", &RangeMeasurementsPython::range_data)
    .add_property("crc", &RangeMeasurementsPython::crc)
  ;

  class_<RawGpsWordPython>("RawGpsWordPython", init<RawGpsWord>())
    .add_property("header", &RawGpsWordPython::header)
    .add_property("prn", &RawGpsWordPython::prn)
    .add_property("nav_word", &RawGpsWordPython::nav_word)
    .add_property("crc", &RawGpsWordPython::crc)
  ;

  class_<GpsEphemerisPython>("GpsEphemerisPython", init<GpsEphemeris>())
    .add_property("header", &GpsEphemerisPython::header)
    .add_property("prn", &GpsEphemerisPython::prn)
    .add_property("time_of_week", &GpsEphemerisPython::time_of_week)
    .add_property("health", &GpsEphemerisPython::health)
    .add_property("issue_of_ephemeris_1",
      &GpsEphemerisPython::issue_of_ephemeris_1)
    .add_property("issue_of_ephemeris_2",
      &GpsEphemerisPython::issue_of_ephemeris_2)
    .add_property("gps_week", &GpsEphemerisPython::gps_week)
    .add_property("z_count_week", &GpsEphemerisPython::z_count_week)
    .add_property("time_of_ephemeris",
      &GpsEphemerisPython::time_of_ephemeris)
    .add_property("semi_major_axis",
      &GpsEphemerisPython::semi_major_axis)
    .add_property("mean_motion_difference",
      &GpsEphemerisPython::mean_motion_difference)
    .add_property("anomoly_reference_time",
      &GpsEphemerisPython::anomoly_reference_time)
    .add_property("eccentricity", &GpsEphemerisPython::eccentricity)
    .add_property("omega", &GpsEphemerisPython::omega)
    .add_property("latitude_cosine",
      &GpsEphemerisPython::latitude_cosine)
    .add_property("latitude_sine", &GpsEphemerisPython::latitude_sine)
    .add_property("orbit_radius_cosine",
      &GpsEphemerisPython::orbit_radius_cosine)
    .add_property("orbit_radius_sine",
      &GpsEphemerisPython::orbit_radius_sine)
    .add_property("inclination_cosine",
      &GpsEphemerisPython::inclination_cosine)
    .add_property("inclination_sine",
      &GpsEphemerisPython::inclination_sine)
    .add_property("inclination_angle",
      &GpsEphemerisPython::inclination_angle)
    .add_property("inclination_angle_rate",
      &GpsEphemerisPython::inclination_angle_rate)
    .add_property("right_ascension",
      &GpsEphemerisPython::right_ascension)
    .add_property("right_ascension_rate",
      &GpsEphemerisPython::right_ascension_rate)
    .add_property("issue_of_data_clock",
      &GpsEphemerisPython::issue_of_data_clock)
    .add_property("sv_clock_correction",
      &GpsEphemerisPython::sv_clock_correction)
    .add_property("group_delay_difference",
      &GpsEphemerisPython::group_delay_difference)
    .add_property("clock_aligning_param_0",
      &GpsEphemerisPython::clock_aligning_param_0)
    .add_property("clock_aligning_param_1",
      &GpsEphemerisPython::clock_aligning_param_1)
    .add_property("clock_aligning_param_2",
      &GpsEphemerisPython::clock_aligning_param_2)
    .add_property("anti_spoofing", &GpsEphemerisPython::anti_spoofing)
    .add_property("corrected_mean_motion",
      &GpsEphemerisPython::corrected_mean_motion)
    .add_property("range_accuracy_variance",
      &GpsEphemerisPython::range_accuracy_variance)
  ;

  class_<AlmanacDataPython>("AlmanacDataPython", init<AlmanacData>())
    .add_property("prn", &AlmanacDataPython::prn)
    .add_property("ref_week", &AlmanacDataPython::ref_week)
    .add_property("ref_time", &AlmanacDataPython::ref_time)
    .add_property("eccentricity", &AlmanacDataPython::eccentricity)
    .add_property("right_ascension_rate",
      &AlmanacDataPython::right_ascension_rate)
    .add_property("right_ascension",
      &AlmanacDataPython::right_ascension)
    .add_property("perigee", &AlmanacDataPython::perigee)
    .add_property("mean_anomoly_of_ref_time",
      &AlmanacDataPython::mean_anomoly_of_ref_time)
    .add_property("clock_aging_param_0",
      &AlmanacDataPython::clock_aging_param_0)
    .add_property("clock_aging_param_1",
      &AlmanacDataPython::clock_aging_param_1)
    .add_property("corrected_mean_motion",
      &AlmanacDataPython::corrected_mean_motion)
    .add_property("semi_major_axis",
      &AlmanacDataPython::semi_major_axis)
    .add_property("inclination_angle",
      &AlmanacDataPython::inclination_angle)
    .add_property("sv_configuration",
      &AlmanacDataPython::sv_configuration)
    .add_property("sv_health", &AlmanacDataPython::sv_health)
    .add_property("sv_health_from_almanac",
      &AlmanacDataPython::sv_health_from_almanac)
    .add_property("anti_spoofing", &AlmanacDataPython::anti_spoofing)
  ;

  class_<AlmanacPython>("AlmanacPython", init<Almanac>())
    .add_property("header", &AlmanacPython::header)
    .add_property("number_of_prns", &AlmanacPython::number_of_prns)
    .add_property("data", &AlmanacPython::data)
    .add_property("crc", &AlmanacPython::crc)
  ;

  class_<BestSatEntry>("BestSatEntry")
    .def_readonly("system", &BestSatEntry::system)
    .def_readonly("satellite_id", &BestSatEntry::satellite_id)
    .def_readonly("status", &BestSatEntry::status)
    .def_readonly("signal_mask", &BestSatEntry::signal_mask)
  ;

  class_<BestSatsPython>("BestSatsPython", init<BestSats>())
    .add_property("header", &BestSatsPython::header)
    .add_property("entries", &BestSatsPython::entries)
    .add_property("sats", &BestSatsPython::sats)
    .add_property("crc", &BestSatsPython::crc)
  ;
}
