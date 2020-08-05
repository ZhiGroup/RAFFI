################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../RaPIDaffin.cpp \
../classifier.cpp \
../dumpable.cpp \
../mapper.cpp \
../ordering.cpp \
../parser.cpp \
../proceed.cpp 

OBJS += \
./RaPIDaffin.o \
./classifier.o \
./dumpable.o \
./mapper.o \
./ordering.o \
./parser.o \
./proceed.o 

CPP_DEPS += \
./RaPIDaffin.d \
./classifier.d \
./dumpable.d \
./mapper.d \
./ordering.d \
./parser.d \
./proceed.d 

CXXFLAGS := -pipe -std=c++14  -Wall  -g

# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ $(CXXFLAGS) -I/home/anaseri/boost/include/ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


