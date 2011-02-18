#ifndef HEMELB_STEERING_BASIC_CLIENTCONNECTION_H
#define HEMELB_STEERING_BASIC_CLIENTCONNECTION_H

namespace hemelb
{
  namespace steering
  {
    class ClientConnection
    {
      public:
        ClientConnection(int iSteeringSessionId);

        int GetWorkingSocket();

        void ReportBroken();

      private:
        static const unsigned int MYPORT = 65250;
        static const unsigned int CONNECTION_BACKLOG = 10;

        int mCurrentSocket;
        bool mIsBroken;
    };
  }
}

#endif /* HEMELB_STEERING_BASIC_CLIENTCONNECTION_H */
