################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../basic_montecarlo.f90 \
../helpers.f90 \
../module_goft.f90 \
../montecarlo.f90 \
../nakajima_zwanzig_shared.f90 \
../numer_fft.f90 \
../numer_fft_row.f90 \
../numer_matrix.f90 \
../std_lapack.f90 \
../std_types.f90 

OBJS += \
./basic_montecarlo.o \
./helpers.o \
./module_goft.o \
./montecarlo.o \
./nakajima_zwanzig_shared.o \
./numer_fft.o \
./numer_fft_row.o \
./numer_matrix.o \
./std_lapack.o \
./std_types.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O3 -Wall -c -fmessage-length=0 -ffree-form -ffree-line-length-none -fopenmp -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

basic_montecarlo.o: ../basic_montecarlo.f90 helpers.o module_goft.o nakajima_zwanzig_shared.o std_types.o

helpers.o: ../helpers.f90 std_types.o

module_goft.o: ../module_goft.f90 helpers.o nakajima_zwanzig_shared.o std_types.o

montecarlo.o: ../montecarlo.f90 basic_montecarlo.o helpers.o std_types.o

nakajima_zwanzig_shared.o: ../nakajima_zwanzig_shared.f90 helpers.o numer_fft.o numer_matrix.o std_lapack.o std_types.o

numer_fft.o: ../numer_fft.f90 std_types.o

numer_fft_row.o: ../numer_fft_row.f90 std_types.o

numer_matrix.o: ../numer_matrix.f90 std_lapack.o std_types.o

std_lapack.o: ../std_lapack.f90 std_types.o

std_types.o: ../std_types.f90


