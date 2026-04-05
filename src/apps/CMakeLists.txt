# src/apps/CMakeLists.txt

add_executable(poisson poisson.cpp)
target_link_libraries(poisson PRIVATE oktal)

# Some tests look for "Poisson"
add_executable(Poisson poisson.cpp)
target_link_libraries(Poisson PRIVATE oktal)

add_executable(tgv tgv.cpp)
target_link_libraries(tgv PRIVATE oktal)
