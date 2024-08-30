#ifndef INTT_ODBC_QUERY_H
#define INTT_ODBC_QUERY_H

class InttOdbcQuery
{
public:
  InttOdbcQuery() = default;
  ~InttOdbcQuery() = default;

  int Verbosity() {return m_verbosity;}
  int Verbosity(int verbosity) {return m_verbosity = verbosity;}

  int Query(int);

  bool IsStreaming();

private:
  static const int m_MAX_NUM_RETRIES = 3000;
  static const int m_MIN_SLEEP_DUR =  200; // milliseconds
  static const int m_MAX_SLEEP_DUR = 3000; // milliseconds

  int m_verbosity{0};
  bool m_query_successful{false};
  bool m_is_streaming{false};
};

#endif//INTT_ODBC_QUERY_H


