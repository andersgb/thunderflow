EXTRA_OPTS = 

PYTHON = /usr/include/python2.7

BOOST_INC = /usr/include
BOOST_LIB = /usr/lib

TARGET = thulip_ext

THULIP_OBJ = api.o load.o core.o thermo_loader.o thermo.o thulip_cpp.o pr_c4.o

$(TARGET).so: $(TARGET).o thulip cpp_interface models
	g++ -shared -Wl,--export-dynamic \
	$(TARGET).o $(THULIP_OBJ) -L$(BOOST_LIB) -lboost_python \
	-L/usr/lib/python2.7/config -lpython2.7 \
	-o $(TARGET).so

$(TARGET).o: $(TARGET).cpp models cpp_interface
	g++ -I$(PYTHON) -I$(BOOST_INC) -Irgrad $(THULIP_OBJ) -ldl -c -fPIC $(TARGET).cpp

all: thunderflow

thulip:
	gcc -c -Wall -fPIC rgrad/thulip/src/*.c

cpp_interface:
	g++ -c -Wall -fPIC -I rgrad -I rgrad/thulip/src rgrad/thulip/cpp/*.cpp

thulip_package_dyn: thulip cpp_interface models
	g++ -shared -Wall -Wl,-soname,libthulippackage.so -o libthulippackage.so $(THULIP_OBJ)

models:
	gcc -c -Wall -fPIC -I rgrad/thulip/src rgrad/models/src/*.c

thunderflow: cpp_interface models
	g++ -Wall -I $(THULIP_OBJ) -ldl main.cpp -o thunderflow $(EXTRA_OPTS)

clean:
	rm -rf *.o thunderflow
