class FileUtils
{
  private:
    // Returns true if the file with the given name exists for reading, false otherwise.
    static bool file_exists(const char * filename);

  public:
    // Exits if the named file doesn't exist or can't be opened for reading.
    static void check_file(const char * filename);
};
