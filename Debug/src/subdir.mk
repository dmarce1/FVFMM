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
../src/node_client.cpp \
../src/node_location.cpp \
../src/node_server.cpp \
../src/node_server_decomp.cpp \
../src/problem.cpp \
../src/roe.cpp \
../src/taylor.cpp 

OBJS += \
./src/grid.o \
./src/grid_fmm.o \
./src/grid_output.o \
./src/lane_emden.o \
./src/main.o \
./src/node_client.o \
./src/node_location.o \
./src/node_server.o \
./src/node_server_decomp.o \
./src/problem.o \
./src/roe.o \
./src/taylor.o 

CPP_DEPS += \
./src/grid.d \
./src/grid_fmm.d \
./src/grid_output.d \
./src/lane_emden.d \
./src/main.d \
./src/node_client.d \
./src/node_location.d \
./src/node_server.d \
./src/node_server_decomp.d \
./src/problem.d \
./src/roe.d \
./src/taylor.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D_GLIBCXX_DEBUG -O0 -g3 -Wall -c -fmessage-length=0  -march=native -std=c++11 `pkg-config --cflags hpx_application_debug` -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


