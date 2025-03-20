## test that changes I make to the rephasing work
add_executable(TestRephase test_rephase.cc)
target_include_directories(TestRephase PUBLIC include)
target_link_libraries(TestRephase PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestPeramSrc test_peram_src.cc)
target_include_directories(TestPeramSrc PUBLIC include)
target_link_libraries(TestPeramSrc PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestNoise test_noise.cc)
target_include_directories(TestNoise PUBLIC include)
target_link_libraries(TestNoise PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestBaryonKernel test_laphBaryonKernel.cc)
target_include_directories(TestBaryonKernel PUBLIC include)
target_link_libraries(TestBaryonKernel PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestBaryonTripletA test_laphBaryonTripletA.cc)
target_include_directories(TestBaryonTripletA PUBLIC include)
target_link_libraries(TestBaryonTripletA PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestBaryonTripletB test_laphBaryonTripletB.cc)
target_include_directories(TestBaryonTripletB PUBLIC include)
target_link_libraries(TestBaryonTripletB PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestCurrentKernel test_laphCurrentKernel.cc)
target_include_directories(TestCurrentKernel PUBLIC include)
target_link_libraries(TestCurrentKernel PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)


