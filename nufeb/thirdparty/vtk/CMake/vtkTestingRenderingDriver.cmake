SET(CMAKE_TESTDRIVER_BEFORE_TESTMAIN
"
    // Set defaults
    vtkTestingInteractor::ValidBaseline = \"Use_-V_for_Baseline\";
    vtkTestingInteractor::TempDirectory =
      std::string(\"${VTK_TEST_OUTPUT_DIR}\");
    vtkTestingInteractor::DataDirectory = std::string(\"Use_-D_for_Data\");

    int interactive = 0;
    for (int ii = 0; ii < ac; ++ii)
      {
      if (strcmp(av[ii], \"-I\") == 0)
        {
        interactive = 1;
        continue;
        }
      if (strcmp(av[ii], \"-V\") == 0 && ii < ac-1)
        {
        vtkTestingInteractor::ValidBaseline = std::string(av[++ii]);
        continue;
        }
      if (strcmp(av[ii], \"-T\") == 0 && ii < ac-1)
        {
        vtkTestingInteractor::TempDirectory = std::string(av[++ii]);
        continue;
        }
      if (strcmp(av[ii], \"-D\") == 0 && ii < ac-1)
        {
        vtkTestingInteractor::DataDirectory = std::string(av[++ii]);
        continue;
        }
      if (strcmp(av[ii], \"-E\") == 0 && ii < ac-1)
        {
        vtkTestingInteractor::ErrorThreshold =
            static_cast<double>(atof(av[++ii]));
        continue;
        }
      }
    vtkSmartPointer<vtkTestingObjectFactory> factory = vtkSmartPointer<vtkTestingObjectFactory>::New();
    if (!interactive)
      {
      // Disable any other overrides before registering our factory.
      vtkObjectFactoryCollection *collection = vtkObjectFactory::GetRegisteredFactories();
      collection->InitTraversal();
      vtkObjectFactory *f = collection->GetNextItem();
      while (f)
        {
        f->Disable(\"vtkRenderWindowInteractor\");
        f = collection->GetNextItem();
        }
      vtkObjectFactory::RegisterFactory(factory);
      }
"
)

SET(CMAKE_TESTDRIVER_AFTER_TESTMAIN
"
   if (result == VTK_SKIP_RETURN_CODE)
     {
     printf(\"Unsupported runtime configuration: Test returned \"
            \"VTK_SKIP_RETURN_CODE. Skipping test.\\n\");
     return result;
     }

   if (!interactive)
     {
     if (vtkTestingInteractor::TestReturnStatus != -1)
        {
        if (vtkTestingInteractor::TestReturnStatus != vtkTesting::PASSED)
          {
          result = EXIT_FAILURE;
          }
        else
          {
          result = EXIT_SUCCESS;
          }
        }
      vtkObjectFactory::UnRegisterFactory(factory);
      }
"
)
