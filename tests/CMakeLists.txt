## test that changes I make to the rephasing work
add_executable(TestRephase test_rephase.cc)
target_include_directories(TestRephase PUBLIC include)
target_link_libraries(TestRephase PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestPeramSrc test_peram_src.cc)
target_include_directories(TestPeramSrc PUBLIC include)
target_link_libraries(TestPeramSrc PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)
