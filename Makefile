PROG =	salsa

SRCS =	driver.f90 driver_input.f90 mo_constants.f90 mo_doctor.f90 \
	mo_exception.f90 mo_kind.f90 mo_salsa.f90 mo_salsa_cloud.f90 \
	mo_salsa_dynamics.f90 mo_salsa_init.f90 mo_salsa_nucleation.f90 \
	mo_salsa_properties.f90 mo_salsa_sizedist.f90 mo_salsa_trac.f90 \
	mo_salsa_update.f90 mo_submctl.f90 mo_time_control.f90

OBJS =	driver.o driver_input.o mo_constants.o mo_doctor.o mo_exception.o \
	mo_kind.o mo_salsa.o mo_salsa_cloud.o mo_salsa_dynamics.o \
	mo_salsa_init.o mo_salsa_nucleation.o mo_salsa_properties.o \
	mo_salsa_sizedist.o mo_salsa_trac.o mo_salsa_update.o mo_submctl.o \
	mo_time_control.o

LIBS =	

CC = cc
CFLAGS = -O
FC = gfortran
FFLAGS = -O
F90 = gfortran
F90FLAGS = -O
LDFLAGS = -s

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

driver.o: driver_input.o mo_kind.o mo_salsa.o mo_salsa_cloud.o \
	mo_salsa_init.o mo_salsa_sizedist.o mo_submctl.o mo_time_control.o
driver_input.o: mo_kind.o
mo_constants.o: mo_kind.o
mo_doctor.o: 
mo_exception.o: mo_doctor.o
mo_kind.o: 
mo_salsa.o: mo_kind.o mo_salsa_dynamics.o mo_salsa_init.o \
	mo_salsa_nucleation.o mo_salsa_properties.o mo_salsa_trac.o \
	mo_salsa_update.o mo_submctl.o
mo_salsa_cloud.o: mo_constants.o mo_kind.o mo_submctl.o
mo_salsa_dynamics.o: mo_constants.o mo_kind.o mo_salsa_init.o \
	mo_salsa_nucleation.o mo_submctl.o mo_time_control.o
mo_salsa_init.o: mo_kind.o mo_submctl.o
mo_salsa_nucleation.o: mo_exception.o mo_kind.o mo_submctl.o
mo_salsa_properties.o: mo_kind.o mo_submctl.o
mo_salsa_sizedist.o: mo_kind.o mo_submctl.o
mo_salsa_trac.o: mo_kind.o mo_submctl.o
mo_salsa_update.o: mo_kind.o mo_submctl.o
mo_submctl.o: mo_kind.o
mo_time_control.o: mo_kind.o
