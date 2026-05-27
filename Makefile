CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -pedantic

TARGET = symplectic_asteroids

SRCS = interface.cc nbody_system.cc orbital_elements.cc units.cc
OBJS = $(SRCS:.cc=.o)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

run: $(TARGET)
	./$(TARGET)