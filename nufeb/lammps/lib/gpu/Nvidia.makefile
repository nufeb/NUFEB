CUDA  = $(NVCC) $(CUDA_INCLUDE) $(CUDA_OPTS) -Icudpp_mini $(CUDA_ARCH) \
             $(CUDA_PRECISION)
CUDR  = $(CUDR_CPP) $(CUDR_OPTS) $(CUDA_PRECISION) $(CUDA_INCLUDE) \
         $(CUDPP_OPT)
CUDA_LINK = $(CUDA_LIB) -lcudart
BIN2C = $(CUDA_HOME)/bin/bin2c

GPU_LIB = $(LIB_DIR)/libgpu.a

# Headers for Geryon
UCL_H  = $(wildcard ./geryon/ucl*.h)
NVC_H  = $(wildcard ./geryon/nvc*.h) $(UCL_H)
NVD_H  = $(wildcard ./geryon/nvd*.h) $(UCL_H) lal_preprocessor.h
# Headers for Pair Stuff
PAIR_H  = lal_atom.h lal_answer.h lal_neighbor_shared.h \
          lal_neighbor.h lal_precision.h lal_device.h \
          lal_balance.h lal_pppm.h

ALL_H = $(NVD_H) $(PAIR_H)

EXECS = $(BIN_DIR)/nvc_get_devices
ifdef CUDPP_OPT
CUDPP = $(OBJ_DIR)/cudpp.o $(OBJ_DIR)/cudpp_plan.o \
        $(OBJ_DIR)/cudpp_maximal_launch.o $(OBJ_DIR)/cudpp_plan_manager.o \
        $(OBJ_DIR)/radixsort_app.cu_o $(OBJ_DIR)/scan_app.cu_o
endif
OBJS = $(OBJ_DIR)/lal_atom.o $(OBJ_DIR)/lal_ans.o \
       $(OBJ_DIR)/lal_neighbor.o $(OBJ_DIR)/lal_neighbor_shared.o \
       $(OBJ_DIR)/lal_device.o $(OBJ_DIR)/lal_base_atomic.o \
       $(OBJ_DIR)/lal_base_charge.o $(OBJ_DIR)/lal_base_ellipsoid.o \
       $(OBJ_DIR)/lal_base_dipole.o $(OBJ_DIR)/lal_base_three.o \
       $(OBJ_DIR)/lal_base_dpd.o \
       $(OBJ_DIR)/lal_pppm.o $(OBJ_DIR)/lal_pppm_ext.o \
       $(OBJ_DIR)/lal_gayberne.o $(OBJ_DIR)/lal_gayberne_ext.o \
       $(OBJ_DIR)/lal_re_squared.o $(OBJ_DIR)/lal_re_squared_ext.o \
       $(OBJ_DIR)/lal_lj.o $(OBJ_DIR)/lal_lj_ext.o \
       $(OBJ_DIR)/lal_lj96.o $(OBJ_DIR)/lal_lj96_ext.o \
       $(OBJ_DIR)/lal_lj_expand.o $(OBJ_DIR)/lal_lj_expand_ext.o \
       $(OBJ_DIR)/lal_lj_coul.o $(OBJ_DIR)/lal_lj_coul_ext.o \
       $(OBJ_DIR)/lal_lj_coul_long.o $(OBJ_DIR)/lal_lj_coul_long_ext.o \
       $(OBJ_DIR)/lal_lj_dsf.o $(OBJ_DIR)/lal_lj_dsf_ext.o \
       $(OBJ_DIR)/lal_lj_class2_long.o $(OBJ_DIR)/lal_lj_class2_long_ext.o \
       $(OBJ_DIR)/lal_coul_long.o $(OBJ_DIR)/lal_coul_long_ext.o \
       $(OBJ_DIR)/lal_morse.o $(OBJ_DIR)/lal_morse_ext.o \
       $(OBJ_DIR)/lal_charmm_long.o $(OBJ_DIR)/lal_charmm_long_ext.o \
       $(OBJ_DIR)/lal_lj_sdk.o $(OBJ_DIR)/lal_lj_sdk_ext.o \
       $(OBJ_DIR)/lal_lj_sdk_long.o $(OBJ_DIR)/lal_lj_sdk_long_ext.o \
       $(OBJ_DIR)/lal_eam.o $(OBJ_DIR)/lal_eam_ext.o \
       $(OBJ_DIR)/lal_eam_fs_ext.o $(OBJ_DIR)/lal_eam_alloy_ext.o \
       $(OBJ_DIR)/lal_buck.o $(OBJ_DIR)/lal_buck_ext.o \
       $(OBJ_DIR)/lal_buck_coul.o $(OBJ_DIR)/lal_buck_coul_ext.o \
       $(OBJ_DIR)/lal_buck_coul_long.o $(OBJ_DIR)/lal_buck_coul_long_ext.o \
       $(OBJ_DIR)/lal_table.o $(OBJ_DIR)/lal_table_ext.o \
       $(OBJ_DIR)/lal_yukawa.o $(OBJ_DIR)/lal_yukawa_ext.o \
       $(OBJ_DIR)/lal_born.o $(OBJ_DIR)/lal_born_ext.o \
       $(OBJ_DIR)/lal_born_coul_wolf.o $(OBJ_DIR)/lal_born_coul_wolf_ext.o \
       $(OBJ_DIR)/lal_born_coul_long.o $(OBJ_DIR)/lal_born_coul_long_ext.o \
       $(OBJ_DIR)/lal_dipole_lj.o $(OBJ_DIR)/lal_dipole_lj_ext.o \
       $(OBJ_DIR)/lal_dipole_lj_sf.o $(OBJ_DIR)/lal_dipole_lj_sf_ext.o \
       $(OBJ_DIR)/lal_colloid.o $(OBJ_DIR)/lal_colloid_ext.o \
       $(OBJ_DIR)/lal_gauss.o $(OBJ_DIR)/lal_gauss_ext.o \
       $(OBJ_DIR)/lal_yukawa_colloid.o $(OBJ_DIR)/lal_yukawa_colloid_ext.o \
       $(OBJ_DIR)/lal_lj_coul_debye.o $(OBJ_DIR)/lal_lj_coul_debye_ext.o \
       $(OBJ_DIR)/lal_coul_dsf.o $(OBJ_DIR)/lal_coul_dsf_ext.o \
       $(OBJ_DIR)/lal_sw.o $(OBJ_DIR)/lal_sw_ext.o \
       $(OBJ_DIR)/lal_vashishta.o $(OBJ_DIR)/lal_vashishta_ext.o \
       $(OBJ_DIR)/lal_beck.o $(OBJ_DIR)/lal_beck_ext.o \
       $(OBJ_DIR)/lal_mie.o $(OBJ_DIR)/lal_mie_ext.o \
       $(OBJ_DIR)/lal_soft.o $(OBJ_DIR)/lal_soft_ext.o \
       $(OBJ_DIR)/lal_lj_coul_msm.o $(OBJ_DIR)/lal_lj_coul_msm_ext.o \
       $(OBJ_DIR)/lal_lj_gromacs.o $(OBJ_DIR)/lal_lj_gromacs_ext.o \
       $(OBJ_DIR)/lal_dpd.o $(OBJ_DIR)/lal_dpd_ext.o \
       $(OBJ_DIR)/lal_tersoff.o $(OBJ_DIR)/lal_tersoff_ext.o \
       $(OBJ_DIR)/lal_tersoff_zbl.o $(OBJ_DIR)/lal_tersoff_zbl_ext.o \
       $(OBJ_DIR)/lal_tersoff_mod.o $(OBJ_DIR)/lal_tersoff_mod_ext.o \
       $(OBJ_DIR)/lal_coul.o $(OBJ_DIR)/lal_coul_ext.o \
       $(OBJ_DIR)/lal_coul_debye.o $(OBJ_DIR)/lal_coul_debye_ext.o \
       $(OBJ_DIR)/lal_zbl.o $(OBJ_DIR)/lal_zbl_ext.o \
       $(OBJ_DIR)/lal_lj_cubic.o $(OBJ_DIR)/lal_lj_cubic_ext.o \
       $(OBJ_DIR)/lal_ufm.o $(OBJ_DIR)/lal_ufm_ext.o \
       $(OBJ_DIR)/lal_dipole_long_lj.o $(OBJ_DIR)/lal_dipole_long_lj_ext.o \
       $(OBJ_DIR)/lal_lj_expand_coul_long.o $(OBJ_DIR)/lal_lj_expand_coul_long_ext.o \
       $(OBJ_DIR)/lal_coul_long_cs.o $(OBJ_DIR)/lal_coul_long_cs_ext.o \
       $(OBJ_DIR)/lal_born_coul_long_cs.o $(OBJ_DIR)/lal_born_coul_long_cs_ext.o \
       $(OBJ_DIR)/lal_born_coul_wolf_cs.o $(OBJ_DIR)/lal_born_coul_wolf_cs_ext.o

CBNS = $(OBJ_DIR)/device.cubin $(OBJ_DIR)/device_cubin.h \
       $(OBJ_DIR)/atom.cubin $(OBJ_DIR)/atom_cubin.h \
       $(OBJ_DIR)/neighbor_cpu.cubin $(OBJ_DIR)/neighbor_cpu_cubin.h \
       $(OBJ_DIR)/neighbor_gpu.cubin $(OBJ_DIR)/neighbor_gpu_cubin.h \
       $(OBJ_DIR)/pppm_f.cubin $(OBJ_DIR)/pppm_f_cubin.h \
       $(OBJ_DIR)/pppm_d.cubin $(OBJ_DIR)/pppm_d_cubin.h \
       $(OBJ_DIR)/ellipsoid_nbor.cubin $(OBJ_DIR)/ellipsoid_nbor_cubin.h \
       $(OBJ_DIR)/gayberne.cubin $(OBJ_DIR)/gayberne_lj.cubin \
       $(OBJ_DIR)/gayberne_cubin.h $(OBJ_DIR)/gayberne_lj_cubin.h \
       $(OBJ_DIR)/re_squared.cubin $(OBJ_DIR)/re_squared_lj.cubin \
       $(OBJ_DIR)/re_squared_cubin.h $(OBJ_DIR)/re_squared_lj_cubin.h \
       $(OBJ_DIR)/lj.cubin $(OBJ_DIR)/lj_cubin.h \
       $(OBJ_DIR)/lj96.cubin $(OBJ_DIR)/lj96_cubin.h \
       $(OBJ_DIR)/lj_expand.cubin $(OBJ_DIR)/lj_expand_cubin.h \
       $(OBJ_DIR)/lj_coul.cubin $(OBJ_DIR)/lj_coul_cubin.h \
       $(OBJ_DIR)/lj_coul_long.cubin $(OBJ_DIR)/lj_coul_long_cubin.h \
       $(OBJ_DIR)/lj_dsf.cubin $(OBJ_DIR)/lj_dsf_cubin.h \
       $(OBJ_DIR)/lj_class2_long.cubin $(OBJ_DIR)/lj_class2_long_cubin.h \
       $(OBJ_DIR)/coul_long.cubin $(OBJ_DIR)/coul_long_cubin.h \
       $(OBJ_DIR)/morse.cubin $(OBJ_DIR)/morse_cubin.h \
       $(OBJ_DIR)/charmm_long.cubin $(OBJ_DIR)/charmm_long_cubin.h \
       $(OBJ_DIR)/lj_sdk.cubin $(OBJ_DIR)/lj_sdk_cubin.h \
       $(OBJ_DIR)/lj_sdk_long.cubin $(OBJ_DIR)/lj_sdk_long_cubin.h \
       $(OBJ_DIR)/eam.cubin $(OBJ_DIR)/eam_cubin.h \
       $(OBJ_DIR)/buck.cubin $(OBJ_DIR)/buck_cubin.h \
       $(OBJ_DIR)/buck_coul_long.cubin $(OBJ_DIR)/buck_coul_long_cubin.h \
       $(OBJ_DIR)/buck_coul.cubin $(OBJ_DIR)/buck_coul_cubin.h \
       $(OBJ_DIR)/table.cubin $(OBJ_DIR)/table_cubin.h \
       $(OBJ_DIR)/yukawa.cubin $(OBJ_DIR)/yukawa_cubin.h \
       $(OBJ_DIR)/born.cubin $(OBJ_DIR)/born_cubin.h \
       $(OBJ_DIR)/born_coul_wolf.cubin $(OBJ_DIR)/born_coul_wolf_cubin.h \
       $(OBJ_DIR)/born_coul_long.cubin $(OBJ_DIR)/born_coul_long_cubin.h \
       $(OBJ_DIR)/dipole_lj.cubin $(OBJ_DIR)/dipole_lj_cubin.h \
       $(OBJ_DIR)/dipole_lj_sf.cubin $(OBJ_DIR)/dipole_lj_sf_cubin.h \
       $(OBJ_DIR)/colloid.cubin $(OBJ_DIR)/colloid_cubin.h \
       $(OBJ_DIR)/gauss.cubin $(OBJ_DIR)/gauss_cubin.h \
       $(OBJ_DIR)/yukawa_colloid.cubin $(OBJ_DIR)/yukawa_colloid_cubin.h \
       $(OBJ_DIR)/lj_coul_debye.cubin $(OBJ_DIR)/lj_coul_debye_cubin.h \
       $(OBJ_DIR)/coul_dsf.cubin $(OBJ_DIR)/coul_dsf_cubin.h \
       $(OBJ_DIR)/sw.cubin $(OBJ_DIR)/sw_cubin.h \
       $(OBJ_DIR)/vashishta.cubin $(OBJ_DIR)/vashishta_cubin.h \
       $(OBJ_DIR)/beck.cubin $(OBJ_DIR)/beck_cubin.h \
       $(OBJ_DIR)/mie.cubin $(OBJ_DIR)/mie_cubin.h \
       $(OBJ_DIR)/soft.cubin $(OBJ_DIR)/soft_cubin.h \
       $(OBJ_DIR)/lj_coul_msm.cubin $(OBJ_DIR)/lj_coul_msm_cubin.h \
       $(OBJ_DIR)/lj_gromacs.cubin $(OBJ_DIR)/lj_gromacs_cubin.h \
       $(OBJ_DIR)/dpd.cubin $(OBJ_DIR)/dpd_cubin.h \
       $(OBJ_DIR)/tersoff.cubin $(OBJ_DIR)/tersoff_cubin.h \
       $(OBJ_DIR)/tersoff_zbl.cubin $(OBJ_DIR)/tersoff_zbl_cubin.h \
       $(OBJ_DIR)/tersoff_mod.cubin $(OBJ_DIR)/tersoff_mod_cubin.h \
       $(OBJ_DIR)/coul.cubin $(OBJ_DIR)/coul_cubin.h \
       $(OBJ_DIR)/coul_debye.cubin $(OBJ_DIR)/coul_debye_cubin.h \
       $(OBJ_DIR)/zbl.cubin $(OBJ_DIR)/zbl_cubin.h \
       $(OBJ_DIR)/lj_cubic.cubin $(OBJ_DIR)/lj_cubic_cubin.h \
       $(OBJ_DIR)/ufm.cubin $(OBJ_DIR)/ufm_cubin.h \
       $(OBJ_DIR)/dipole_long_lj.cubin $(OBJ_DIR)/dipole_long_lj_cubin.h \
       $(OBJ_DIR)/lj_expand_coul_long.cubin $(OBJ_DIR)/lj_expand_coul_long_cubin.h \
       $(OBJ_DIR)/coul_long_cs.cubin $(OBJ_DIR)/coul_long_cs_cubin.h \
       $(OBJ_DIR)/born_coul_long_cs.cubin $(OBJ_DIR)/born_coul_long_cs_cubin.h \
       $(OBJ_DIR)/born_coul_wolf_cs.cubin $(OBJ_DIR)/born_coul_wolf_cs_cubin.h

all: $(OBJ_DIR) $(GPU_LIB) $(EXECS)

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/cudpp.o: cudpp_mini/cudpp.cpp
	$(CUDR) -o $@ -c cudpp_mini/cudpp.cpp -Icudpp_mini

$(OBJ_DIR)/cudpp_plan.o: cudpp_mini/cudpp_plan.cpp
	$(CUDR) -o $@ -c cudpp_mini/cudpp_plan.cpp -Icudpp_mini

$(OBJ_DIR)/cudpp_maximal_launch.o: cudpp_mini/cudpp_maximal_launch.cpp
	$(CUDR) -o $@ -c cudpp_mini/cudpp_maximal_launch.cpp -Icudpp_mini

$(OBJ_DIR)/cudpp_plan_manager.o: cudpp_mini/cudpp_plan_manager.cpp
	$(CUDR) -o $@ -c cudpp_mini/cudpp_plan_manager.cpp -Icudpp_mini

$(OBJ_DIR)/radixsort_app.cu_o: cudpp_mini/radixsort_app.cu
	$(CUDA) -o $@ -c cudpp_mini/radixsort_app.cu

$(OBJ_DIR)/scan_app.cu_o: cudpp_mini/scan_app.cu
	$(CUDA) -o $@ -c cudpp_mini/scan_app.cu

$(OBJ_DIR)/atom.cubin: lal_atom.cu lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_atom.cu

$(OBJ_DIR)/atom_cubin.h: $(OBJ_DIR)/atom.cubin
	$(BIN2C) -c -n atom $(OBJ_DIR)/atom.cubin > $(OBJ_DIR)/atom_cubin.h

$(OBJ_DIR)/lal_atom.o: lal_atom.cpp lal_atom.h $(NVD_H) $(OBJ_DIR)/atom_cubin.h
	$(CUDR) -o $@ -c lal_atom.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_ans.o: lal_answer.cpp lal_answer.h $(NVD_H)
	$(CUDR) -o $@ -c lal_answer.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/neighbor_cpu.cubin: lal_neighbor_cpu.cu lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_neighbor_cpu.cu

$(OBJ_DIR)/neighbor_cpu_cubin.h: $(OBJ_DIR)/neighbor_cpu.cubin
	$(BIN2C) -c -n neighbor_cpu $(OBJ_DIR)/neighbor_cpu.cubin > $(OBJ_DIR)/neighbor_cpu_cubin.h

$(OBJ_DIR)/neighbor_gpu.cubin: lal_neighbor_gpu.cu lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_neighbor_gpu.cu

$(OBJ_DIR)/neighbor_gpu_cubin.h: $(OBJ_DIR)/neighbor_gpu.cubin
	$(BIN2C) -c -n neighbor_gpu $(OBJ_DIR)/neighbor_gpu.cubin > $(OBJ_DIR)/neighbor_gpu_cubin.h

$(OBJ_DIR)/lal_neighbor_shared.o: lal_neighbor_shared.cpp lal_neighbor_shared.h $(OBJ_DIR)/neighbor_cpu_cubin.h $(OBJ_DIR)/neighbor_gpu_cubin.h $(NVD_H)
	$(CUDR) -o $@ -c lal_neighbor_shared.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_neighbor.o: lal_neighbor.cpp lal_neighbor.h lal_neighbor_shared.h $(NVD_H)
	$(CUDR) -o $@ -c lal_neighbor.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/device.cubin: lal_device.cu lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_device.cu

$(OBJ_DIR)/device_cubin.h: $(OBJ_DIR)/device.cubin
	$(BIN2C) -c -n device $(OBJ_DIR)/device.cubin > $(OBJ_DIR)/device_cubin.h

$(OBJ_DIR)/lal_device.o: lal_device.cpp lal_device.h $(ALL_H) $(OBJ_DIR)/device_cubin.h
	$(CUDR) -o $@ -c lal_device.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_base_atomic.o: $(ALL_H) lal_base_atomic.h lal_base_atomic.cpp
	$(CUDR) -o $@ -c lal_base_atomic.cpp

$(OBJ_DIR)/lal_base_charge.o: $(ALL_H) lal_base_charge.h lal_base_charge.cpp
	$(CUDR) -o $@ -c lal_base_charge.cpp

$(OBJ_DIR)/lal_base_ellipsoid.o: $(ALL_H) lal_base_ellipsoid.h lal_base_ellipsoid.cpp $(OBJ_DIR)/ellipsoid_nbor_cubin.h
	$(CUDR) -o $@ -c lal_base_ellipsoid.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_base_dipole.o: $(ALL_H) lal_base_dipole.h lal_base_dipole.cpp
	$(CUDR) -o $@ -c lal_base_dipole.cpp

$(OBJ_DIR)/lal_base_three.o: $(ALL_H) lal_base_three.h lal_base_three.cpp
	$(CUDR) -o $@ -c lal_base_three.cpp

$(OBJ_DIR)/lal_base_dpd.o: $(ALL_H) lal_base_dpd.h lal_base_dpd.cpp
	$(CUDR) -o $@ -c lal_base_dpd.cpp

$(OBJ_DIR)/pppm_f.cubin: lal_pppm.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -Dgrdtyp=float -Dgrdtyp4=float4 -o $@ lal_pppm.cu

$(OBJ_DIR)/pppm_f_cubin.h: $(OBJ_DIR)/pppm_f.cubin
	$(BIN2C) -c -n pppm_f $(OBJ_DIR)/pppm_f.cubin > $(OBJ_DIR)/pppm_f_cubin.h

$(OBJ_DIR)/pppm_d.cubin: lal_pppm.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -Dgrdtyp=double -Dgrdtyp4=double4 -o $@ lal_pppm.cu

$(OBJ_DIR)/pppm_d_cubin.h: $(OBJ_DIR)/pppm_d.cubin
	$(BIN2C) -c -n pppm_d $(OBJ_DIR)/pppm_d.cubin > $(OBJ_DIR)/pppm_d_cubin.h

$(OBJ_DIR)/lal_pppm.o: $(ALL_H) lal_pppm.h lal_pppm.cpp $(OBJ_DIR)/pppm_f_cubin.h $(OBJ_DIR)/pppm_d_cubin.h
	$(CUDR) -o $@ -c lal_pppm.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_pppm_ext.o: $(ALL_H) lal_pppm.h lal_pppm_ext.cpp
	$(CUDR) -o $@ -c lal_pppm_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/ellipsoid_nbor.cubin: lal_ellipsoid_nbor.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_ellipsoid_nbor.cu

$(OBJ_DIR)/ellipsoid_nbor_cubin.h: $(OBJ_DIR)/ellipsoid_nbor.cubin
	$(BIN2C) -c -n ellipsoid_nbor $(OBJ_DIR)/ellipsoid_nbor.cubin > $(OBJ_DIR)/ellipsoid_nbor_cubin.h

$(OBJ_DIR)/gayberne.cubin: lal_gayberne.cu lal_precision.h lal_ellipsoid_extra.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_gayberne.cu

$(OBJ_DIR)/gayberne_lj.cubin: lal_gayberne_lj.cu lal_precision.h lal_ellipsoid_extra.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_gayberne_lj.cu

$(OBJ_DIR)/gayberne_cubin.h: $(OBJ_DIR)/gayberne.cubin
	$(BIN2C) -c -n gayberne $(OBJ_DIR)/gayberne.cubin > $(OBJ_DIR)/gayberne_cubin.h

$(OBJ_DIR)/gayberne_lj_cubin.h: $(OBJ_DIR)/gayberne_lj.cubin
	$(BIN2C) -c -n gayberne_lj $(OBJ_DIR)/gayberne_lj.cubin > $(OBJ_DIR)/gayberne_lj_cubin.h

$(OBJ_DIR)/lal_gayberne.o: $(ALL_H) lal_gayberne.h lal_gayberne.cpp $(OBJ_DIR)/gayberne_cubin.h $(OBJ_DIR)/gayberne_lj_cubin.h $(OBJ_DIR)/lal_base_ellipsoid.o
	$(CUDR) -o $@ -c lal_gayberne.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_gayberne_ext.o: $(ALL_H) $(OBJ_DIR)/lal_gayberne.o lal_gayberne_ext.cpp
	$(CUDR) -o $@ -c lal_gayberne_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/re_squared.cubin: lal_re_squared.cu lal_precision.h lal_ellipsoid_extra.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_re_squared.cu

$(OBJ_DIR)/re_squared_lj.cubin: lal_re_squared_lj.cu lal_precision.h lal_ellipsoid_extra.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_re_squared_lj.cu

$(OBJ_DIR)/re_squared_cubin.h: $(OBJ_DIR)/re_squared.cubin
	$(BIN2C) -c -n re_squared $(OBJ_DIR)/re_squared.cubin > $(OBJ_DIR)/re_squared_cubin.h

$(OBJ_DIR)/re_squared_lj_cubin.h: $(OBJ_DIR)/re_squared_lj.cubin
	$(BIN2C) -c -n re_squared_lj $(OBJ_DIR)/re_squared_lj.cubin > $(OBJ_DIR)/re_squared_lj_cubin.h

$(OBJ_DIR)/lal_re_squared.o: $(ALL_H) lal_re_squared.h lal_re_squared.cpp $(OBJ_DIR)/re_squared_cubin.h $(OBJ_DIR)/re_squared_lj_cubin.h $(OBJ_DIR)/lal_base_ellipsoid.o
	$(CUDR) -o $@ -c lal_re_squared.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_re_squared_ext.o: $(ALL_H) $(OBJ_DIR)/lal_re_squared.o lal_re_squared_ext.cpp
	$(CUDR) -o $@ -c lal_re_squared_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj.cubin: lal_lj.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_lj.cu

$(OBJ_DIR)/lj_cubin.h: $(OBJ_DIR)/lj.cubin $(OBJ_DIR)/lj.cubin
	$(BIN2C) -c -n lj $(OBJ_DIR)/lj.cubin > $(OBJ_DIR)/lj_cubin.h

$(OBJ_DIR)/lal_lj.o: $(ALL_H) lal_lj.h lal_lj.cpp $(OBJ_DIR)/lj_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_lj.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_ext.o: $(ALL_H) lal_lj.h lal_lj_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_lj_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul.cubin: lal_lj_coul.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_lj_coul.cu

$(OBJ_DIR)/lj_coul_cubin.h: $(OBJ_DIR)/lj_coul.cubin $(OBJ_DIR)/lj_coul.cubin
	$(BIN2C) -c -n lj_coul $(OBJ_DIR)/lj_coul.cubin > $(OBJ_DIR)/lj_coul_cubin.h

$(OBJ_DIR)/lal_lj_coul.o: $(ALL_H) lal_lj_coul.h lal_lj_coul.cpp $(OBJ_DIR)/lj_coul_cubin.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_lj_coul.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_coul_ext.o: $(ALL_H) lal_lj_coul.h lal_lj_coul_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_lj_coul_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_class2_long.cubin: lal_lj_class2_long.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_lj_class2_long.cu

$(OBJ_DIR)/lj_class2_long_cubin.h: $(OBJ_DIR)/lj_class2_long.cubin $(OBJ_DIR)/lj_class2_long.cubin
	$(BIN2C) -c -n lj_class2_long $(OBJ_DIR)/lj_class2_long.cubin > $(OBJ_DIR)/lj_class2_long_cubin.h

$(OBJ_DIR)/lal_lj_class2_long.o: $(ALL_H) lal_lj_class2_long.h lal_lj_class2_long.cpp $(OBJ_DIR)/lj_class2_long_cubin.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_lj_class2_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_class2_long_ext.o: $(ALL_H) lal_lj_class2_long.h lal_lj_class2_long_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_lj_class2_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/coul_long.cubin: lal_coul_long.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_coul_long.cu

$(OBJ_DIR)/coul_long_cubin.h: $(OBJ_DIR)/coul_long.cubin $(OBJ_DIR)/coul_long.cubin
	$(BIN2C) -c -n coul_long $(OBJ_DIR)/coul_long.cubin > $(OBJ_DIR)/coul_long_cubin.h

$(OBJ_DIR)/lal_coul_long.o: $(ALL_H) lal_coul_long.h lal_coul_long.cpp $(OBJ_DIR)/coul_long_cubin.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_coul_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_coul_long_ext.o: $(ALL_H) lal_coul_long.h lal_coul_long_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_coul_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul_long.cubin: lal_lj_coul_long.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_lj_coul_long.cu

$(OBJ_DIR)/lj_coul_long_cubin.h: $(OBJ_DIR)/lj_coul_long.cubin $(OBJ_DIR)/lj_coul_long.cubin
	$(BIN2C) -c -n lj_coul_long $(OBJ_DIR)/lj_coul_long.cubin > $(OBJ_DIR)/lj_coul_long_cubin.h

$(OBJ_DIR)/lal_lj_coul_long.o: $(ALL_H) lal_lj_coul_long.h lal_lj_coul_long.cpp $(OBJ_DIR)/lj_coul_long_cubin.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_lj_coul_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_coul_long_ext.o: $(ALL_H) lal_lj_coul_long.h lal_lj_coul_long_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_lj_coul_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_dsf.cubin: lal_lj_dsf.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_lj_dsf.cu

$(OBJ_DIR)/lj_dsf_cubin.h: $(OBJ_DIR)/lj_dsf.cubin $(OBJ_DIR)/lj_dsf.cubin
	$(BIN2C) -c -n lj_dsf $(OBJ_DIR)/lj_dsf.cubin > $(OBJ_DIR)/lj_dsf_cubin.h

$(OBJ_DIR)/lal_lj_dsf.o: $(ALL_H) lal_lj_dsf.h lal_lj_dsf.cpp $(OBJ_DIR)/lj_dsf_cubin.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_lj_dsf.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_dsf_ext.o: $(ALL_H) lal_lj_dsf.h lal_lj_dsf_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_lj_dsf_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/morse.cubin: lal_morse.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_morse.cu

$(OBJ_DIR)/morse_cubin.h: $(OBJ_DIR)/morse.cubin $(OBJ_DIR)/morse.cubin
	$(BIN2C) -c -n morse $(OBJ_DIR)/morse.cubin > $(OBJ_DIR)/morse_cubin.h

$(OBJ_DIR)/lal_morse.o: $(ALL_H) lal_morse.h lal_morse.cpp $(OBJ_DIR)/morse_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_morse.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_morse_ext.o: $(ALL_H) lal_morse.h lal_morse_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_morse_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/charmm_long.cubin: lal_charmm_long.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_charmm_long.cu

$(OBJ_DIR)/charmm_long_cubin.h: $(OBJ_DIR)/charmm_long.cubin $(OBJ_DIR)/charmm_long.cubin
	$(BIN2C) -c -n charmm_long $(OBJ_DIR)/charmm_long.cubin > $(OBJ_DIR)/charmm_long_cubin.h

$(OBJ_DIR)/lal_charmm_long.o: $(ALL_H) lal_charmm_long.h lal_charmm_long.cpp $(OBJ_DIR)/charmm_long_cubin.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_charmm_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_charmm_long_ext.o: $(ALL_H) lal_charmm_long.h lal_charmm_long_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_charmm_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj96.cubin: lal_lj96.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_lj96.cu

$(OBJ_DIR)/lj96_cubin.h: $(OBJ_DIR)/lj96.cubin $(OBJ_DIR)/lj96.cubin
	$(BIN2C) -c -n lj96 $(OBJ_DIR)/lj96.cubin > $(OBJ_DIR)/lj96_cubin.h

$(OBJ_DIR)/lal_lj96.o: $(ALL_H) lal_lj96.h lal_lj96.cpp $(OBJ_DIR)/lj96_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_lj96.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj96_ext.o: $(ALL_H) lal_lj96.h lal_lj96_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_lj96_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_expand.cubin: lal_lj_expand.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_lj_expand.cu

$(OBJ_DIR)/lj_expand_cubin.h: $(OBJ_DIR)/lj_expand.cubin $(OBJ_DIR)/lj_expand.cubin
	$(BIN2C) -c -n lj_expand $(OBJ_DIR)/lj_expand.cubin > $(OBJ_DIR)/lj_expand_cubin.h

$(OBJ_DIR)/lal_lj_expand.o: $(ALL_H) lal_lj_expand.h lal_lj_expand.cpp $(OBJ_DIR)/lj_expand_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_lj_expand.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_expand_ext.o: $(ALL_H) lal_lj_expand.h lal_lj_expand_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_lj_expand_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_sdk.cubin: lal_lj_sdk.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_lj_sdk.cu

$(OBJ_DIR)/lj_sdk_cubin.h: $(OBJ_DIR)/lj_sdk.cubin $(OBJ_DIR)/lj_sdk.cubin
	$(BIN2C) -c -n lj_sdk $(OBJ_DIR)/lj_sdk.cubin > $(OBJ_DIR)/lj_sdk_cubin.h

$(OBJ_DIR)/lal_lj_sdk.o: $(ALL_H) lal_lj_sdk.h lal_lj_sdk.cpp $(OBJ_DIR)/lj_sdk_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_lj_sdk.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_sdk_ext.o: $(ALL_H) lal_lj_sdk.h lal_lj_sdk_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_lj_sdk_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_sdk_long.cubin: lal_lj_sdk_long.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_lj_sdk_long.cu

$(OBJ_DIR)/lj_sdk_long_cubin.h: $(OBJ_DIR)/lj_sdk_long.cubin $(OBJ_DIR)/lj_sdk_long.cubin
	$(BIN2C) -c -n lj_sdk_long $(OBJ_DIR)/lj_sdk_long.cubin > $(OBJ_DIR)/lj_sdk_long_cubin.h

$(OBJ_DIR)/lal_lj_sdk_long.o: $(ALL_H) lal_lj_sdk_long.h lal_lj_sdk_long.cpp $(OBJ_DIR)/lj_sdk_long_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_lj_sdk_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_sdk_long_ext.o: $(ALL_H) lal_lj_sdk_long.h lal_lj_sdk_long_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_lj_sdk_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/eam.cubin: lal_eam.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_eam.cu

$(OBJ_DIR)/eam_cubin.h: $(OBJ_DIR)/eam.cubin $(OBJ_DIR)/eam.cubin
	$(BIN2C) -c -n eam $(OBJ_DIR)/eam.cubin > $(OBJ_DIR)/eam_cubin.h

$(OBJ_DIR)/lal_eam.o: $(ALL_H) lal_eam.h lal_eam.cpp $(OBJ_DIR)/eam_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_eam.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_eam_ext.o: $(ALL_H) lal_eam.h lal_eam_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_eam_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_eam_fs_ext.o: $(ALL_H) lal_eam.h lal_eam_fs_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_eam_fs_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_eam_alloy_ext.o: $(ALL_H) lal_eam.h lal_eam_alloy_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_eam_alloy_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/buck.cubin: lal_buck.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_buck.cu

$(OBJ_DIR)/buck_cubin.h: $(OBJ_DIR)/buck.cubin $(OBJ_DIR)/buck.cubin
	$(BIN2C) -c -n buck $(OBJ_DIR)/buck.cubin > $(OBJ_DIR)/buck_cubin.h

$(OBJ_DIR)/lal_buck.o: $(ALL_H) lal_buck.h lal_buck.cpp $(OBJ_DIR)/buck_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_buck.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_buck_ext.o: $(ALL_H) lal_buck.h lal_buck_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_buck_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/buck_coul.cubin: lal_buck_coul.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_buck_coul.cu

$(OBJ_DIR)/buck_coul_cubin.h: $(OBJ_DIR)/buck_coul.cubin $(OBJ_DIR)/buck_coul.cubin
	$(BIN2C) -c -n buck_coul $(OBJ_DIR)/buck_coul.cubin > $(OBJ_DIR)/buck_coul_cubin.h

$(OBJ_DIR)/lal_buck_coul.o: $(ALL_H) lal_buck_coul.h lal_buck_coul.cpp $(OBJ_DIR)/buck_coul_cubin.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_buck_coul.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_buck_coul_ext.o: $(ALL_H) lal_buck_coul.h lal_buck_coul_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_buck_coul_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/buck_coul_long.cubin: lal_buck_coul_long.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_buck_coul_long.cu

$(OBJ_DIR)/buck_coul_long_cubin.h: $(OBJ_DIR)/buck_coul_long.cubin $(OBJ_DIR)/buck_coul_long.cubin
	$(BIN2C) -c -n buck_coul_long $(OBJ_DIR)/buck_coul_long.cubin > $(OBJ_DIR)/buck_coul_long_cubin.h

$(OBJ_DIR)/lal_buck_coul_long.o: $(ALL_H) lal_buck_coul_long.h lal_buck_coul_long.cpp $(OBJ_DIR)/buck_coul_long_cubin.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_buck_coul_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_buck_coul_long_ext.o: $(ALL_H) lal_buck_coul_long.h lal_buck_coul_long_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_buck_coul_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/table.cubin: lal_table.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_table.cu

$(OBJ_DIR)/table_cubin.h: $(OBJ_DIR)/table.cubin $(OBJ_DIR)/table.cubin
	$(BIN2C) -c -n table $(OBJ_DIR)/table.cubin > $(OBJ_DIR)/table_cubin.h

$(OBJ_DIR)/lal_table.o: $(ALL_H) lal_table.h lal_table.cpp $(OBJ_DIR)/table_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_table.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_table_ext.o: $(ALL_H) lal_table.h lal_table_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_table_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/yukawa.cubin: lal_yukawa.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_yukawa.cu

$(OBJ_DIR)/yukawa_cubin.h: $(OBJ_DIR)/yukawa.cubin $(OBJ_DIR)/yukawa.cubin
	$(BIN2C) -c -n yukawa $(OBJ_DIR)/yukawa.cubin > $(OBJ_DIR)/yukawa_cubin.h

$(OBJ_DIR)/lal_yukawa.o: $(ALL_H) lal_yukawa.h lal_yukawa.cpp $(OBJ_DIR)/yukawa_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_yukawa.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_yukawa_ext.o: $(ALL_H) lal_yukawa.h lal_yukawa_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_yukawa_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/born.cubin: lal_born.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_born.cu

$(OBJ_DIR)/born_cubin.h: $(OBJ_DIR)/born.cubin $(OBJ_DIR)/born.cubin
	$(BIN2C) -c -n born $(OBJ_DIR)/born.cubin > $(OBJ_DIR)/born_cubin.h

$(OBJ_DIR)/lal_born.o: $(ALL_H) lal_born.h lal_born.cpp $(OBJ_DIR)/born_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_born.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_born_ext.o: $(ALL_H) lal_born.h lal_born_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_born_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/born_coul_wolf.cubin: lal_born_coul_wolf.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_born_coul_wolf.cu

$(OBJ_DIR)/born_coul_wolf_cubin.h: $(OBJ_DIR)/born_coul_wolf.cubin $(OBJ_DIR)/born_coul_wolf.cubin
	$(BIN2C) -c -n born_coul_wolf $(OBJ_DIR)/born_coul_wolf.cubin > $(OBJ_DIR)/born_coul_wolf_cubin.h

$(OBJ_DIR)/lal_born_coul_wolf.o: $(ALL_H) lal_born_coul_wolf.h lal_born_coul_wolf.cpp $(OBJ_DIR)/born_coul_wolf_cubin.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_born_coul_wolf.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_born_coul_wolf_ext.o: $(ALL_H) lal_born_coul_wolf.h lal_born_coul_wolf_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_born_coul_wolf_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/born_coul_long.cubin: lal_born_coul_long.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_born_coul_long.cu

$(OBJ_DIR)/born_coul_long_cubin.h: $(OBJ_DIR)/born_coul_long.cubin $(OBJ_DIR)/born_coul_long.cubin
	$(BIN2C) -c -n born_coul_long $(OBJ_DIR)/born_coul_long.cubin > $(OBJ_DIR)/born_coul_long_cubin.h

$(OBJ_DIR)/lal_born_coul_long.o: $(ALL_H) lal_born_coul_long.h lal_born_coul_long.cpp $(OBJ_DIR)/born_coul_long_cubin.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_born_coul_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_born_coul_long_ext.o: $(ALL_H) lal_born_coul_long.h lal_born_coul_long_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_born_coul_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/dipole_lj.cubin: lal_dipole_lj.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_dipole_lj.cu

$(OBJ_DIR)/dipole_lj_cubin.h: $(OBJ_DIR)/dipole_lj.cubin $(OBJ_DIR)/dipole_lj.cubin
	$(BIN2C) -c -n dipole_lj $(OBJ_DIR)/dipole_lj.cubin > $(OBJ_DIR)/dipole_lj_cubin.h

$(OBJ_DIR)/lal_dipole_lj.o: $(ALL_H) lal_dipole_lj.h lal_dipole_lj.cpp $(OBJ_DIR)/dipole_lj_cubin.h $(OBJ_DIR)/lal_base_dipole.o
	$(CUDR) -o $@ -c lal_dipole_lj.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_dipole_lj_ext.o: $(ALL_H) lal_dipole_lj.h lal_dipole_lj_ext.cpp lal_base_dipole.h
	$(CUDR) -o $@ -c lal_dipole_lj_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/dipole_lj_sf.cubin: lal_dipole_lj_sf.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_dipole_lj_sf.cu

$(OBJ_DIR)/dipole_lj_sf_cubin.h: $(OBJ_DIR)/dipole_lj_sf.cubin $(OBJ_DIR)/dipole_lj_sf.cubin
	$(BIN2C) -c -n dipole_lj_sf $(OBJ_DIR)/dipole_lj_sf.cubin > $(OBJ_DIR)/dipole_lj_sf_cubin.h

$(OBJ_DIR)/lal_dipole_lj_sf.o: $(ALL_H) lal_dipole_lj_sf.h lal_dipole_lj_sf.cpp $(OBJ_DIR)/dipole_lj_sf_cubin.h $(OBJ_DIR)/lal_base_dipole.o
	$(CUDR) -o $@ -c lal_dipole_lj_sf.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_dipole_lj_sf_ext.o: $(ALL_H) lal_dipole_lj_sf.h lal_dipole_lj_sf_ext.cpp lal_base_dipole.h
	$(CUDR) -o $@ -c lal_dipole_lj_sf_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/colloid.cubin: lal_colloid.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_colloid.cu

$(OBJ_DIR)/colloid_cubin.h: $(OBJ_DIR)/colloid.cubin $(OBJ_DIR)/colloid.cubin
	$(BIN2C) -c -n colloid $(OBJ_DIR)/colloid.cubin > $(OBJ_DIR)/colloid_cubin.h

$(OBJ_DIR)/lal_colloid.o: $(ALL_H) lal_colloid.h lal_colloid.cpp $(OBJ_DIR)/colloid_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_colloid.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_colloid_ext.o: $(ALL_H) lal_colloid.h lal_colloid_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_colloid_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/gauss.cubin: lal_gauss.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_gauss.cu

$(OBJ_DIR)/gauss_cubin.h: $(OBJ_DIR)/gauss.cubin $(OBJ_DIR)/gauss.cubin
	$(BIN2C) -c -n gauss $(OBJ_DIR)/gauss.cubin > $(OBJ_DIR)/gauss_cubin.h

$(OBJ_DIR)/lal_gauss.o: $(ALL_H) lal_gauss.h lal_gauss.cpp $(OBJ_DIR)/gauss_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_gauss.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_gauss_ext.o: $(ALL_H) lal_gauss.h lal_gauss_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_gauss_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/yukawa_colloid.cubin: lal_yukawa_colloid.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_yukawa_colloid.cu

$(OBJ_DIR)/yukawa_colloid_cubin.h: $(OBJ_DIR)/yukawa_colloid.cubin $(OBJ_DIR)/yukawa_colloid.cubin
	$(BIN2C) -c -n yukawa_colloid $(OBJ_DIR)/yukawa_colloid.cubin > $(OBJ_DIR)/yukawa_colloid_cubin.h

$(OBJ_DIR)/lal_yukawa_colloid.o: $(ALL_H) lal_yukawa_colloid.h lal_yukawa_colloid.cpp $(OBJ_DIR)/yukawa_colloid_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_yukawa_colloid.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_yukawa_colloid_ext.o: $(ALL_H) lal_yukawa_colloid.h lal_yukawa_colloid_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_yukawa_colloid_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul_debye.cubin: lal_lj_coul_debye.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_lj_coul_debye.cu

$(OBJ_DIR)/lj_coul_debye_cubin.h: $(OBJ_DIR)/lj_coul_debye.cubin $(OBJ_DIR)/lj_coul_debye.cubin
	$(BIN2C) -c -n lj_coul_debye $(OBJ_DIR)/lj_coul_debye.cubin > $(OBJ_DIR)/lj_coul_debye_cubin.h

$(OBJ_DIR)/lal_lj_coul_debye.o: $(ALL_H) lal_lj_coul_debye.h lal_lj_coul_debye.cpp $(OBJ_DIR)/lj_coul_debye_cubin.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_lj_coul_debye.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_coul_debye_ext.o: $(ALL_H) lal_lj_coul_debye.h lal_lj_coul_debye_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_lj_coul_debye_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/coul_dsf.cubin: lal_coul_dsf.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_coul_dsf.cu

$(OBJ_DIR)/coul_dsf_cubin.h: $(OBJ_DIR)/coul_dsf.cubin $(OBJ_DIR)/coul_dsf.cubin
	$(BIN2C) -c -n coul_dsf $(OBJ_DIR)/coul_dsf.cubin > $(OBJ_DIR)/coul_dsf_cubin.h

$(OBJ_DIR)/lal_coul_dsf.o: $(ALL_H) lal_coul_dsf.h lal_coul_dsf.cpp $(OBJ_DIR)/coul_dsf_cubin.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_coul_dsf.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_coul_dsf_ext.o: $(ALL_H) lal_coul_dsf.h lal_coul_dsf_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_coul_dsf_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/vashishta.cubin: lal_vashishta.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_vashishta.cu

$(OBJ_DIR)/vashishta_cubin.h: $(OBJ_DIR)/vashishta.cubin $(OBJ_DIR)/vashishta.cubin
	$(BIN2C) -c -n vashishta $(OBJ_DIR)/vashishta.cubin > $(OBJ_DIR)/vashishta_cubin.h

$(OBJ_DIR)/lal_vashishta.o: $(ALL_H) lal_vashishta.h lal_vashishta.cpp $(OBJ_DIR)/vashishta_cubin.h $(OBJ_DIR)/lal_base_three.o
	$(CUDR) -o $@ -c lal_vashishta.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_vashishta_ext.o: $(ALL_H) lal_vashishta.h lal_vashishta_ext.cpp lal_base_three.h
	$(CUDR) -o $@ -c lal_vashishta_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/sw.cubin: lal_sw.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_sw.cu

$(OBJ_DIR)/sw_cubin.h: $(OBJ_DIR)/sw.cubin $(OBJ_DIR)/sw.cubin
	$(BIN2C) -c -n sw $(OBJ_DIR)/sw.cubin > $(OBJ_DIR)/sw_cubin.h

$(OBJ_DIR)/lal_sw.o: $(ALL_H) lal_sw.h lal_sw.cpp $(OBJ_DIR)/sw_cubin.h $(OBJ_DIR)/lal_base_three.o
	$(CUDR) -o $@ -c lal_sw.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_sw_ext.o: $(ALL_H) lal_sw.h lal_sw_ext.cpp lal_base_three.h
	$(CUDR) -o $@ -c lal_sw_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/beck.cubin: lal_beck.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_beck.cu

$(OBJ_DIR)/beck_cubin.h: $(OBJ_DIR)/beck.cubin $(OBJ_DIR)/beck.cubin
	$(BIN2C) -c -n beck $(OBJ_DIR)/beck.cubin > $(OBJ_DIR)/beck_cubin.h

$(OBJ_DIR)/lal_beck.o: $(ALL_H) lal_beck.h lal_beck.cpp $(OBJ_DIR)/beck_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_beck.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_beck_ext.o: $(ALL_H) lal_beck.h lal_beck_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_beck_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/mie.cubin: lal_mie.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_mie.cu

$(OBJ_DIR)/mie_cubin.h: $(OBJ_DIR)/mie.cubin $(OBJ_DIR)/mie.cubin
	$(BIN2C) -c -n mie $(OBJ_DIR)/mie.cubin > $(OBJ_DIR)/mie_cubin.h

$(OBJ_DIR)/lal_mie.o: $(ALL_H) lal_mie.h lal_mie.cpp $(OBJ_DIR)/mie_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_mie.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_mie_ext.o: $(ALL_H) lal_mie.h lal_mie_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_mie_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/soft.cubin: lal_soft.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_soft.cu

$(OBJ_DIR)/soft_cubin.h: $(OBJ_DIR)/soft.cubin $(OBJ_DIR)/soft.cubin
	$(BIN2C) -c -n soft $(OBJ_DIR)/soft.cubin > $(OBJ_DIR)/soft_cubin.h

$(OBJ_DIR)/lal_soft.o: $(ALL_H) lal_soft.h lal_soft.cpp $(OBJ_DIR)/soft_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_soft.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_soft_ext.o: $(ALL_H) lal_soft.h lal_soft_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_soft_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_coul_msm.cubin: lal_lj_coul_msm.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_lj_coul_msm.cu

$(OBJ_DIR)/lj_coul_msm_cubin.h: $(OBJ_DIR)/lj_coul_msm.cubin $(OBJ_DIR)/lj_coul_msm.cubin
	$(BIN2C) -c -n lj_coul_msm $(OBJ_DIR)/lj_coul_msm.cubin > $(OBJ_DIR)/lj_coul_msm_cubin.h

$(OBJ_DIR)/lal_lj_coul_msm.o: $(ALL_H) lal_lj_coul_msm.h lal_lj_coul_msm.cpp $(OBJ_DIR)/lj_coul_msm_cubin.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_lj_coul_msm.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_coul_msm_ext.o: $(ALL_H) lal_lj_coul_msm.h lal_lj_coul_msm_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_lj_coul_msm_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_gromacs.cubin: lal_lj_gromacs.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_lj_gromacs.cu

$(OBJ_DIR)/lj_gromacs_cubin.h: $(OBJ_DIR)/lj_gromacs.cubin $(OBJ_DIR)/lj_gromacs.cubin
	$(BIN2C) -c -n lj_gromacs $(OBJ_DIR)/lj_gromacs.cubin > $(OBJ_DIR)/lj_gromacs_cubin.h

$(OBJ_DIR)/lal_lj_gromacs.o: $(ALL_H) lal_lj_gromacs.h lal_lj_gromacs.cpp $(OBJ_DIR)/lj_gromacs_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_lj_gromacs.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_gromacs_ext.o: $(ALL_H) lal_lj_gromacs.h lal_lj_gromacs_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_lj_gromacs_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/dpd.cubin: lal_dpd.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_dpd.cu

$(OBJ_DIR)/dpd_cubin.h: $(OBJ_DIR)/dpd.cubin $(OBJ_DIR)/dpd.cubin
	$(BIN2C) -c -n dpd $(OBJ_DIR)/dpd.cubin > $(OBJ_DIR)/dpd_cubin.h

$(OBJ_DIR)/ufm.cubin: lal_ufm.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_ufm.cu

$(OBJ_DIR)/ufm_cubin.h: $(OBJ_DIR)/ufm.cubin $(OBJ_DIR)/ufm.cubin
	$(BIN2C) -c -n ufm $(OBJ_DIR)/ufm.cubin > $(OBJ_DIR)/ufm_cubin.h

$(OBJ_DIR)/lal_ufm.o: $(ALL_H) lal_ufm.h lal_ufm.cpp $(OBJ_DIR)/ufm_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_ufm.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_ufm_ext.o: $(ALL_H) lal_ufm.h lal_ufm_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_ufm_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_dpd.o: $(ALL_H) lal_dpd.h lal_dpd.cpp $(OBJ_DIR)/dpd_cubin.h $(OBJ_DIR)/lal_base_dpd.o
	$(CUDR) -o $@ -c lal_dpd.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_dpd_ext.o: $(ALL_H) lal_dpd.h lal_dpd_ext.cpp lal_base_dpd.h
	$(CUDR) -o $@ -c lal_dpd_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/tersoff.cubin: lal_tersoff.cu lal_precision.h lal_tersoff_extra.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_tersoff.cu

$(OBJ_DIR)/tersoff_cubin.h: $(OBJ_DIR)/tersoff.cubin $(OBJ_DIR)/tersoff.cubin
	$(BIN2C) -c -n tersoff $(OBJ_DIR)/tersoff.cubin > $(OBJ_DIR)/tersoff_cubin.h

$(OBJ_DIR)/lal_tersoff.o: $(ALL_H) lal_tersoff.h lal_tersoff.cpp $(OBJ_DIR)/tersoff_cubin.h $(OBJ_DIR)/lal_base_three.o
	$(CUDR) -o $@ -c lal_tersoff.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_tersoff_ext.o: $(ALL_H) lal_tersoff.h lal_tersoff_ext.cpp lal_base_three.h
	$(CUDR) -o $@ -c lal_tersoff_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/tersoff_zbl.cubin: lal_tersoff_zbl.cu lal_precision.h lal_tersoff_zbl_extra.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_tersoff_zbl.cu

$(OBJ_DIR)/tersoff_zbl_cubin.h: $(OBJ_DIR)/tersoff_zbl.cubin $(OBJ_DIR)/tersoff_zbl.cubin
	$(BIN2C) -c -n tersoff_zbl $(OBJ_DIR)/tersoff_zbl.cubin > $(OBJ_DIR)/tersoff_zbl_cubin.h

$(OBJ_DIR)/lal_tersoff_zbl.o: $(ALL_H) lal_tersoff_zbl.h lal_tersoff_zbl.cpp $(OBJ_DIR)/tersoff_zbl_cubin.h $(OBJ_DIR)/lal_base_three.o
	$(CUDR) -o $@ -c lal_tersoff_zbl.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_tersoff_zbl_ext.o: $(ALL_H) lal_tersoff_zbl.h lal_tersoff_zbl_ext.cpp lal_base_three.h
	$(CUDR) -o $@ -c lal_tersoff_zbl_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/tersoff_mod.cubin: lal_tersoff_mod.cu lal_precision.h lal_tersoff_mod_extra.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_tersoff_mod.cu

$(OBJ_DIR)/tersoff_mod_cubin.h: $(OBJ_DIR)/tersoff_mod.cubin $(OBJ_DIR)/tersoff_mod.cubin
	$(BIN2C) -c -n tersoff_mod $(OBJ_DIR)/tersoff_mod.cubin > $(OBJ_DIR)/tersoff_mod_cubin.h

$(OBJ_DIR)/lal_tersoff_mod.o: $(ALL_H) lal_tersoff_mod.h lal_tersoff_mod.cpp $(OBJ_DIR)/tersoff_mod_cubin.h $(OBJ_DIR)/lal_base_three.o
	$(CUDR) -o $@ -c lal_tersoff_mod.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_tersoff_mod_ext.o: $(ALL_H) lal_tersoff_mod.h lal_tersoff_mod_ext.cpp lal_base_three.h
	$(CUDR) -o $@ -c lal_tersoff_mod_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/coul.cubin: lal_coul.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_coul.cu

$(OBJ_DIR)/coul_cubin.h: $(OBJ_DIR)/coul.cubin $(OBJ_DIR)/coul.cubin
	$(BIN2C) -c -n coul $(OBJ_DIR)/coul.cubin > $(OBJ_DIR)/coul_cubin.h

$(OBJ_DIR)/lal_coul.o: $(ALL_H) lal_coul.h lal_coul.cpp $(OBJ_DIR)/coul_cubin.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_coul.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_coul_ext.o: $(ALL_H) lal_coul.h lal_coul_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_coul_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/coul_debye.cubin: lal_coul_debye.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_coul_debye.cu

$(OBJ_DIR)/coul_debye_cubin.h: $(OBJ_DIR)/coul_debye.cubin $(OBJ_DIR)/coul_debye.cubin
	$(BIN2C) -c -n coul_debye $(OBJ_DIR)/coul_debye.cubin > $(OBJ_DIR)/coul_debye_cubin.h

$(OBJ_DIR)/lal_coul_debye.o: $(ALL_H) lal_coul_debye.h lal_coul_debye.cpp $(OBJ_DIR)/coul_debye_cubin.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_coul_debye.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_coul_debye_ext.o: $(ALL_H) lal_coul_debye.h lal_coul_debye_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_coul_debye_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/zbl.cubin: lal_zbl.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_zbl.cu

$(OBJ_DIR)/zbl_cubin.h: $(OBJ_DIR)/zbl.cubin $(OBJ_DIR)/zbl.cubin
	$(BIN2C) -c -n zbl $(OBJ_DIR)/zbl.cubin > $(OBJ_DIR)/zbl_cubin.h

$(OBJ_DIR)/lal_zbl.o: $(ALL_H) lal_zbl.h lal_zbl.cpp $(OBJ_DIR)/zbl_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_zbl.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_zbl_ext.o: $(ALL_H) lal_zbl.h lal_zbl_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_zbl_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_cubic.cubin: lal_lj_cubic.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_lj_cubic.cu

$(OBJ_DIR)/lj_cubic_cubin.h: $(OBJ_DIR)/lj_cubic.cubin $(OBJ_DIR)/lj_cubic.cubin
	$(BIN2C) -c -n lj_cubic $(OBJ_DIR)/lj_cubic.cubin > $(OBJ_DIR)/lj_cubic_cubin.h

$(OBJ_DIR)/lal_lj_cubic.o: $(ALL_H) lal_lj_cubic.h lal_lj_cubic.cpp $(OBJ_DIR)/lj_cubic_cubin.h $(OBJ_DIR)/lal_base_atomic.o
	$(CUDR) -o $@ -c lal_lj_cubic.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_cubic_ext.o: $(ALL_H) lal_lj_cubic.h lal_lj_cubic_ext.cpp lal_base_atomic.h
	$(CUDR) -o $@ -c lal_lj_cubic_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/dipole_long_lj.cubin: lal_dipole_long_lj.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_dipole_long_lj.cu

$(OBJ_DIR)/dipole_long_lj_cubin.h: $(OBJ_DIR)/dipole_long_lj.cubin $(OBJ_DIR)/dipole_long_lj.cubin
	$(BIN2C) -c -n dipole_long_lj $(OBJ_DIR)/dipole_long_lj.cubin > $(OBJ_DIR)/dipole_long_lj_cubin.h

$(OBJ_DIR)/lal_dipole_long_lj.o: $(ALL_H) lal_dipole_long_lj.h lal_dipole_long_lj.cpp $(OBJ_DIR)/dipole_long_lj_cubin.h $(OBJ_DIR)/lal_base_dipole.o
	$(CUDR) -o $@ -c lal_dipole_long_lj.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_dipole_long_lj_ext.o: $(ALL_H) lal_dipole_long_lj.h lal_dipole_long_lj_ext.cpp lal_base_dipole.h
	$(CUDR) -o $@ -c lal_dipole_long_lj_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lj_expand_coul_long.cubin: lal_lj_expand_coul_long.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_lj_expand_coul_long.cu

$(OBJ_DIR)/lj_expand_coul_long_cubin.h: $(OBJ_DIR)/lj_expand_coul_long.cubin $(OBJ_DIR)/lj_expand_coul_long.cubin
	$(BIN2C) -c -n lj_expand_coul_long $(OBJ_DIR)/lj_expand_coul_long.cubin > $(OBJ_DIR)/lj_expand_coul_long_cubin.h

$(OBJ_DIR)/lal_lj_expand_coul_long.o: $(ALL_H) lal_lj_expand_coul_long.h lal_lj_expand_coul_long.cpp $(OBJ_DIR)/lj_expand_coul_long_cubin.h $(OBJ_DIR)/lal_base_charge.o
	$(CUDR) -o $@ -c lal_lj_expand_coul_long.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_lj_expand_coul_long_ext.o: $(ALL_H) lal_lj_expand_coul_long.h lal_lj_expand_coul_long_ext.cpp lal_base_charge.h
	$(CUDR) -o $@ -c lal_lj_expand_coul_long_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/coul_long_cs.cubin: lal_coul_long_cs.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_coul_long_cs.cu

$(OBJ_DIR)/coul_long_cs_cubin.h: $(OBJ_DIR)/coul_long_cs.cubin $(OBJ_DIR)/coul_long_cs.cubin
	$(BIN2C) -c -n coul_long_cs $(OBJ_DIR)/coul_long_cs.cubin > $(OBJ_DIR)/coul_long_cs_cubin.h

$(OBJ_DIR)/lal_coul_long_cs.o: $(ALL_H) lal_coul_long_cs.h lal_coul_long_cs.cpp $(OBJ_DIR)/coul_long_cs_cubin.h $(OBJ_DIR)/lal_base_charge.o $(OBJ_DIR)/lal_coul_long.o
	$(CUDR) -o $@ -c lal_coul_long_cs.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_coul_long_cs_ext.o: $(ALL_H) lal_coul_long_cs.h lal_coul_long_cs_ext.cpp lal_coul_long.h
	$(CUDR) -o $@ -c lal_coul_long_cs_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/born_coul_long_cs.cubin: lal_born_coul_long_cs.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_born_coul_long_cs.cu

$(OBJ_DIR)/born_coul_long_cs_cubin.h: $(OBJ_DIR)/born_coul_long_cs.cubin $(OBJ_DIR)/born_coul_long_cs.cubin
	$(BIN2C) -c -n born_coul_long_cs $(OBJ_DIR)/born_coul_long_cs.cubin > $(OBJ_DIR)/born_coul_long_cs_cubin.h

$(OBJ_DIR)/lal_born_coul_long_cs.o: $(ALL_H) lal_born_coul_long_cs.h lal_born_coul_long_cs.cpp $(OBJ_DIR)/born_coul_long_cs_cubin.h $(OBJ_DIR)/lal_base_charge.o $(OBJ_DIR)/lal_born_coul_long.o
	$(CUDR) -o $@ -c lal_born_coul_long_cs.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_born_coul_long_cs_ext.o: $(ALL_H) lal_born_coul_long_cs.h lal_born_coul_long_cs_ext.cpp lal_born_coul_long.h
	$(CUDR) -o $@ -c lal_born_coul_long_cs_ext.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/born_coul_wolf_cs.cubin: lal_born_coul_wolf_cs.cu lal_precision.h lal_preprocessor.h
	$(CUDA) --cubin -DNV_KERNEL -o $@ lal_born_coul_wolf_cs.cu

$(OBJ_DIR)/born_coul_wolf_cs_cubin.h: $(OBJ_DIR)/born_coul_wolf_cs.cubin $(OBJ_DIR)/born_coul_wolf_cs.cubin
	$(BIN2C) -c -n born_coul_wolf_cs $(OBJ_DIR)/born_coul_wolf_cs.cubin > $(OBJ_DIR)/born_coul_wolf_cs_cubin.h

$(OBJ_DIR)/lal_born_coul_wolf_cs.o: $(ALL_H) lal_born_coul_wolf_cs.h lal_born_coul_wolf_cs.cpp $(OBJ_DIR)/born_coul_wolf_cs_cubin.h $(OBJ_DIR)/lal_base_charge.o $(OBJ_DIR)/lal_born_coul_wolf.o
	$(CUDR) -o $@ -c lal_born_coul_wolf_cs.cpp -I$(OBJ_DIR)

$(OBJ_DIR)/lal_born_coul_wolf_cs_ext.o: $(ALL_H) lal_born_coul_wolf_cs.h lal_born_coul_wolf_cs_ext.cpp lal_born_coul_wolf.h
	$(CUDR) -o $@ -c lal_born_coul_wolf_cs_ext.cpp -I$(OBJ_DIR)

$(BIN_DIR)/nvc_get_devices: ./geryon/ucl_get_devices.cpp $(NVD_H)
	$(CUDR) -o $@ ./geryon/ucl_get_devices.cpp -DUCL_CUDADR $(CUDA_LIB) -lcuda

$(GPU_LIB): $(OBJS) $(CUDPP)
	$(AR) -crusv $(GPU_LIB) $(OBJS) $(CUDPP)
	@cp $(EXTRAMAKE) Makefile.lammps

clean:
	-rm -f $(EXECS) $(GPU_LIB) $(OBJS) $(CUDPP) $(CBNS) *.linkinfo

veryclean: clean
	-rm -rf *~ *.linkinfo

cleanlib:
	-rm -f $(EXECS) $(GPU_LIB) $(OBJS) $(CBNS) *.linkinfo
