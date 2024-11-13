
include_directories(${QudaLaph}/include/)

target_sources(quda_laph PRIVATE
  quda_laph.cc
  byte_handler.cc
  dilution_handler.cc
  dilution_scheme_info.cc
  field_ops.cc
  field_smearing_info.cc
  filelist_info.cc
  gauge_configuration_handler.cc
  gauge_configuration_info.cc
  gluon_smearing_handler.cc
  inverter_info.cc
  io_fm_handler.cc
  laph_eigen_info.cc
  laph_noise.cc
  laph_noise_info.cc
  laph_stdio.cc
  latt_field.cc
  layout_info.cc
  perambulator_handler.cc
  quark_action_info.cc
  quark_handler.cc
  quark_smearing_handler.cc
  quda_info.cc
  read_gauge_field.cc
  stop_watch.cc
  task_quark_line_ends.cc
  task_quark_perambulators.cc
  task_smear_gauge_field.cc
  task_smear_quark_field.cc
  utils.cc
  xml_handler.cc
)
