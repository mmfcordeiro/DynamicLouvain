cmake_minimum_required (VERSION 2.8)

project (community-latest)

add_executable (convert main_convert.cpp graph.h graph.cpp) 
add_executable (community main_community.cpp graph_binary.h graph_binary.cpp community.h community.cpp) 
add_executable (community_sandbox main_community_sandbox.cpp graph_binary.h graph_binary.cpp community.h community.cpp)
add_executable (hierarchy main_hierarchy.cpp)
add_executable (random main_random.cpp) 
 