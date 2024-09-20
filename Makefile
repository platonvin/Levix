.ONESHELL:

#setting up include and lib directories for dependencies
I = -Isrc
L = 

lib_objs := \

ifeq ($(OS),Windows_NT)
	REQUIRED_LIBS =
	STATIC_OR_DYNAMIC += -static
else
	REQUIRED_LIBS =
endif
	
always_enabled_flags = -fno-exceptions -Wuninitialized -std=c++20
special_otp_flags = -O3

srcs := \
	src/main.cpp\

profile = -D_PRINTLINE -DVKNDEBUG
profile = 
#default target
all: init main

obj/main.o: src/main.cpp init vcpkg_installed_eval
	c++ $(special_otp_flags) $(always_enabled_flags) $(I) $(args) $(profile) -MMD -MP -c $< -o $@ 
DEPS = $(obj/main.d)
-include $(DEPS)

main: init vcpkg_installed_eval library build_main
ifeq ($(OS),Windows_NT)
	.\main
else
	./main
endif

library: init vcpkg_installed_eval $(lib_objs)
	ar rvs lib/liblevix.a $(lib_objs)

build_main: library obj/main.o
	c++ -o main obj/main.o -Llib $(flags) $(I) $(L) $(REQUIRED_LIBS) $(STATIC_OR_DYNAMIC)

init: obj lib
obj:
ifeq ($(OS),Windows_NT)
	mkdir "obj"
else
	mkdir -p obj
endif

lib:
ifeq ($(OS),Windows_NT)
	mkdir "lib"
else
	mkdir -p lib
endif

clean:
ifeq ($(OS),Windows_NT)
	del "obj\*.o" "lib\*.a" "examples\*.spv"
else
	rm -rf obj/*.o lib/*.a examples/*.spv
endif

.PHONY: vcpkg_installed_eval 
vcpkg_installed_eval: vcpkg_installed
	$(eval OTHER_DIRS := $(filter-out vcpkg_installed/vcpkg, $(wildcard vcpkg_installed/*)) )
	$(eval INCLUDE_LIST := $(addsuffix /include, $(OTHER_DIRS)) )
	$(eval LIB_LIST := $(addsuffix /lib, $(OTHER_DIRS)) )
	$(eval I += $(addprefix -I, $(INCLUDE_LIST)) )
	$(eval L += $(addprefix -L, $(LIB_LIST)) )

vcpkg_installed:
	echo installind vcpkg dependencies. Please do not interrupt
	vcpkg install