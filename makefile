DEBUG = -g
EXTERN_LIB = -L./lib -lgmp -lcrypto -lsodium -ldl -lpthread
INCLUDE = -I./include/ 
CXX = g++
TARGET1 = dragonfly_sta
TARGET2 = dragonfly_ap
SRC1 = $(wildcard ./src/*.cpp)
SRC2 = dragonfly_sta.cpp
SRC3 = dragonfly_ap.cpp
OBJ1 = $(patsubst %.cpp, %.o, $(SRC1))
OBJ2 = $(patsubst %.cpp, %.o, $(SRC2))
OBJ3 = $(patsubst %.cpp, %.o, $(SRC3))
SRC = ./alg/

all: $(TARGET1) $(TARGET2)

$(TARGET1):$(OBJ1) $(OBJ2)
	$(CXX) $^ -o $@ $(EXTERN_LIB)

$(TARGET2):$(OBJ1) $(OBJ3)
	$(CXX) $^ -o $@ $(EXTERN_LIB)

#编译SRC变量代表的目录下的.cpp文件
%.o:$(SRC)%.cpp
	$(CXX) $(DEBUG) -std=c++11 -c $< -o $@ $(INCLUDE)

#编译当前目录下的.cpp文件
%.o:%.cpp
	$(CXX) $(DEBUG) -std=c++11 -c $< -o $@ $(INCLUDE)

#防止外面有clean文件，阻止执行clean
.PHONY:clean

#防止外面有clean文件，阻止执行clean
.PHONY:clean

clean:
	-rm -rf $(TARGET) $(OBJ1) $(OBJ2) $(OBJ3) $(TARGET1) $(TARGET2)
