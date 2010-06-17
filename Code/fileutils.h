class FileUtils
{
  private:
    // Returns true if the file with the given name exists for reading, false otherwise.
    static bool file_exists(const char * filename);

    // Function to select directory contents that are not "." or ".."
    // The hack is necessary because of different versions of 'scandir'
    // used by different compilers.
    #ifdef DARWIN
      static int SelectOnlyContents (struct direct *entry);
    #else
      static int SelectOnlyContents (const struct direct *entry);
    #endif

  public:
    // Exits if the named file doesn't exist or can't be opened for reading.
    static void check_file(const char * filename);

    // Delete all files within a directory.
    static int DeleteDirContents (char *pathname);

    // Function to create the directory of given path, which user group and anyone
    // can read write and execute.
    static void MakeDirAllRXW(char* dirPath);
};
