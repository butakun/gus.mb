// $Id: Clock.h 41 2010-07-06 10:48:07Z kato $
#ifndef INCLUDED_CLOCK_H__
#define INCLUDED_CLOCK_H__

class Clock
{
public:
    static void Initialize(double t0 = 0.0);
    static Clock* GetInstance();

    double Time() const;
    void Advance(double dt);

protected:
    Clock(double t0 = 0.0);

private:
    static Clock* mInstance;

    double mT;
};

#endif // INCLUDED_CLOCK_H__

