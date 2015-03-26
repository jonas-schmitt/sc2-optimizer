// This header was taken from the NuSiF exercise skeleton!
// Why rewrite :)

#ifndef DEBUG_H
#define DEBUG_H

#include <sstream>

//===================================================================================================================
//
//  CHECK Macro used for tests, is activated in Debug and Release Mode
//
//===================================================================================================================

#define CHECK_MSG(X, MSG) \
   if( !(X) )  \
   { std::stringstream ss; \
     ss << MSG; \
     internal::checkFct ( (X), #X, ss.str(), __FILE__, __LINE__ );\
   }

#define CHECK(X) \
   if( !(X) ) {  internal::checkFct ( (X), #X, "", __FILE__, __LINE__ ); }




#define WARN(MSG) \
   { std::stringstream ss; \
     ss << MSG; \
     internal::warnFct ( ss.str(), __FILE__, __LINE__ );\
   }




//===================================================================================================================
//
//  ASSERT Macro checks the given expression in Debug mode, disabled in Release mode
//  PROGRESS Macro prints out feedback on the order of execution of the program in Debug mode, disabled in Release mode
//
//===================================================================================================================


#ifndef NDEBUG

#define ASSERT_MSG(X, MSG) \
   if( !(X) )  \
   { std::stringstream ss; \
     ss << MSG; \
     internal::assertFct ( (X), #X, ss.str(), __FILE__, __LINE__ );\
   }

#define ASSERT(X) \
   if( !(X) ) {  internal::assertFct ( (X), #X, "", __FILE__, __LINE__ ); }


/*#define PROGRESS(MSG) \
    { \
        std::stringstream ss;  \
        ss << MSG; \
        internal::progressFct(ss.str()); \
    }
*/
#define PROGRESS(MSG)

#else

#define ASSERT_MSG(X, MSG)
#define ASSERT(X)
#define PROGRESS(MSG)

#endif //NDEBUG







namespace internal
{

   void checkFct ( bool b, const char * const expression, const std::string & message,
                   const char * const filename, int line );
   void assertFct( bool b, const char * const expression, const std::string & message,
                   const char * const filename, int line );

   void warnFct( const std::string & message,
                 const char * const filename, int line );

   void progressFct(const std::string & message);

}



#endif // DEBUG_H
