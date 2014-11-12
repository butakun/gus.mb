// $Id: Profiler.h 189 2012-01-19 07:40:08Z kato $
#ifndef INCLUDED_PROFILER_H__
#define INCLUDED_PROFILER_H__

class Profiler
{
public:
    static Profiler* GetInstance();

    virtual ~Profiler();

    double CheckPoint(); // returns the elapsed time since the last checkpoint.

protected:
    Profiler();

private:
    static Profiler* mInstance;

    double mTimeStamp;
};

#endif // INCLUDED_PROFILER_H__

