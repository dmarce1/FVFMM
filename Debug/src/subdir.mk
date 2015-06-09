################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/grid.cpp \
../src/grid_fmm.cpp \
../src/grid_output.cpp \
../src/lane_emden.cpp \
../src/main.cpp \
../src/problem.cpp \
../src/roe.cpp \
../src/taylor.cpp 

OBJS += \
./src/grid.o \
./src/grid_fmm.o \
./src/grid_output.o \
./src/lane_emden.o \
./src/main.o \
./src/problem.o \
./src/roe.o \
./src/taylor.o 

CPP_DEPS += \
./src/grid.d \
./src/grid_fmm.d \
./src/grid_output.d \
./src/lane_emden.d \
./src/main.d \
./src/problem.d \
./src/roe.d \
./src/taylor.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++11 -D_GLIBCXX_DEBUG -I/home/dmarce1/include -I/usr/include/mpi -O0 -g3 -Wall -c -fmessage-length=0 `pkg-config --cflags hpx_application_debug`   -march=native -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


