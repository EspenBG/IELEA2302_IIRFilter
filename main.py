# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig

show_figures = False  # Display the figure when the program is run
save_figures = True  # Save the figures to the same folder as the program
print_coefficient = True  # Enable the printing of coefficients for the filters


def plot_bode(omega, amp, fig_num=None, title=""):
    plt.close(fig_num)
    bode = plt.figure(fig_num)
    bode.suptitle(title)

    # Make the subplot for the amplitude
    plot_amplitude = bode.add_subplot(211)
    plot_amplitude.semilogx(omega, 20 * np.log10(abs(amp)))
    #plot_amplitude.set_xlim(0, 50)
    #plot_amplitude.set_ylim(-100, 10)
    # plot_amplitude.xlabel('Frequency [radians / second]')
    plot_amplitude.set_ylabel('Amplitude [dB]')
    plot_amplitude.margins(0, 0.1)
    plot_amplitude.grid(which='both', axis='both')
    #plot_amplitude.axvline(100, color='green')  # cutoff frequency

    # Make the subplot for the phase response
    phase = bode.add_subplot(212, sharex=plot_amplitude)
    phase.semilogx(omega, np.angle(amp))
    # phase.title('Butterworth filter frequency response')
    phase.set_xlabel('Frequency [radians / second]')
    phase.set_ylabel('Phase [radians]')
    phase.margins(0, 0.1)
    phase.grid(which='both', axis='both')
    #phase.axvline(100, color='green')  # cutoff frequency
    bode.show()


def compare2disc(analog_filter, ts=0.1, title="", fig_num=[None, None], step_time=[0, 10]):
    #ts = 1/100

    # Make filters in the s-plane


    w, h = sig.freqs(analog_filter[0], analog_filter[1], worN=np.logspace(-2, 3, num=10000))

    # Plot the bode and step responses for the filters
    #plot_bode(w, h, 1, "e∂")

    # Transform the filters to the z-domain
    # Using zoh
    filter_z_zoh = sig.cont2discrete(system=analog_filter, dt=ts, method="zoh") # bilinear, zoh
    filter_z_tustin = sig.cont2discrete(system=analog_filter, dt=ts, method="bilinear") # bilinear, zoh
    filter_z_warp = sig.cont2discrete(system=analog_filter, dt=ts, method="zoh") # bilinear, zoh
    # using tustin

    if print_coefficient:
        print(title, " analog: \n b:", analog_filter[0], "\n a: ", analog_filter[1])
        print(title, " discrete zoh: \n b:", filter_z_zoh[0][0], "\n a: ", filter_z_zoh[1])
        print(title, " discrete tustin: \n b:", filter_z_tustin[0][0], "\n a: ", filter_z_tustin[1])

    # using tustin with warping

    # Plot the z-domain representation of the filter (bode and step respons)


    w_z_zoh, h_z_zoh = sig.freqz(filter_z_zoh[0][0], filter_z_zoh[1], fs=1/filter_z_zoh[2], worN=10000)
    w_z_tustin, h_z_tustin = sig.freqz(filter_z_tustin[0][0], filter_z_tustin[1], fs=1/filter_z_tustin[2], worN=10000)

    #w_z_zoh, h_z_zoh = sig.freqz([0.00995], [1, -0.99], fs=100, worN=5000)
    #plt.close(fig_num)
    plot_bode = plt.figure(dpi=200, figsize=(12.8, 7), num=fig_num[0])
    plot_bode.suptitle("Bodeplot "+ title)

    # Make the subplot for the amplitude
    plot_amplitude = plot_bode.add_subplot(211)
    plot_amplitude.semilogx(w, 20 * np.log10(abs(h)), alpha=1, label=r'$\vert H(s) \vert$ Continuous', linewidth=2)
    plot_amplitude.semilogx(w_z_tustin * 2 * np.pi, 20 * np.log10(abs(h_z_tustin)), alpha=1, label=r"$\vert H(z) \vert$ Discrete tustin", linewidth=1)
    plot_amplitude.semilogx(w_z_zoh*2*np.pi, 20 * np.log10(abs(h_z_zoh)), "g", alpha=1, label=r"$\vert H(z) \vert$ Discrete zoh", linewidth=1)
    plot_amplitude.axvline(1/ts*np.pi, color='red', ls=":")  # cutoff frequency
    plot_amplitude.legend()
    plot_amplitude.set_xlim(0.01, 400)
    plot_amplitude.set_ylim(-80, 5)
    # plot_amplitude.xlabel('Frequency [radians / second]')
    plot_amplitude.set_ylabel('Amplitude [dB]')
    plot_amplitude.margins(0, 0.1)
    plot_amplitude.grid(which='both', axis='both')
    #plot_amplitude.axvline(100, color='green')  # cutoff frequency

    # Make the subplot for the phase response
    plot_phase = plot_bode.add_subplot(212, sharex=plot_amplitude)
    plot_phase.semilogx(w, np.angle(h), label=r"$\angle H(s)$ Continuous", linewidth=2)
    plot_phase.semilogx(w_z_tustin*2*np.pi, np.unwrap(np.angle(h_z_tustin)), label=r"$\angle H(z)$ Discrete tustin", linewidth=1)
    plot_phase.semilogx(w_z_zoh*2*np.pi, np.unwrap(np.angle(h_z_zoh)), label=r"$\angle H(z)$ Discrete zoh", linewidth=1)
    plot_phase.axvline(1/ts*np.pi, color='red', ls=":")  # cutoff frequency
    plot_phase.legend()
    # plot_phase.set_ylim(-np.pi, np.pi)
    # phase.title('Butterworth filter frequency response')
    plot_phase.set_xlabel('Frequency [radians / second]')
    plot_phase.set_ylabel('Phase [radians]')
    plot_phase.margins(0, 0.1)
    plot_phase.grid(which='both', axis='both')
    #phase.axvline(100, color='green')  # cutoff frequency
    #bode.show()
    plot_bode.tight_layout()

#    plot_bode(w_z_zoh, h_z_zoh
    #print()

    # make the stepresponses
    step_length = int((step_time[1] - step_time[0])/ts)
    t = np.linspace(step_time[0], step_time[1], step_length, endpoint=True)
    # Stepsignal
    u_n = np.ones(step_length)
    # Make the steprespns for discrete zoh system
    s_n_zoh = sig.lfilter(filter_z_zoh[0][0], filter_z_zoh[1], u_n)
    s_n_tustin = sig.lfilter(filter_z_tustin[0][0], filter_z_tustin[1], u_n)
    s_t_cont = sig.lsim(analog_filter, u_n, t)

    # make the impulse response
    # Make the impulse for discrete zoh system
    h_n_zoh = sig.dimpulse(filter_z_zoh, t=t)
    h_n_tustin = sig.dimpulse(filter_z_tustin, t=t)
    h_t_cont = sig.impulse(analog_filter, T=t)

    plot_step = plt.figure(dpi=200, figsize=(12.8, 7), num=fig_num[1])
    plot_step.suptitle("Step and impulse response for " + title)

    # Make the subplot for the stepresponse
    ax1_step = plot_step.add_subplot(211)
    ax1_step.set_title("Step response")
    ax1_step.plot(s_t_cont[0], s_t_cont[1], "C0", alpha=1, label=r'$ U(s)\cdot H(s) $ Continuous', linewidth=2)
    ax1_step.plot(t, s_n_tustin, "C1", alpha=1, label=r'$U(z)\cdot H(z) $ Discrete tustin', linewidth=1)
    ax1_step.plot(t, s_n_zoh, "C2", alpha=1, label=r'$ U(z)\cdot H(z)$ Discrete zoh', linewidth=1)

    #ax1_step.axvline(1 / ts * np.pi, color='red', ls=":")  # cutoff frequency
    ax1_step.legend()
    #ax1_step.set_xlim(0.01, 400)
    #ax1_step.set_ylim(-80, 5)
    ax1_step.set_xlabel('Time [seconds]')
    ax1_step.set_ylabel('Amplitude')
    ax1_step.margins(0, 0.1)
    ax1_step.grid(which='both', axis='both')


    # Make the subplot for the impulse response
    ax1_impulse = plot_step.add_subplot(212, sharex=ax1_step)
    ax1_impulse.set_title("Impulse response")
    ax1_impulse.plot(h_t_cont[0], h_t_cont[1], "C0", alpha=1, label=r'$ U(s)\cdot H(s) $ Continuous', linewidth=2)
    ax1_impulse.plot(h_n_tustin[0], h_n_tustin[1][0], "C1", alpha=1, label=r'$U(z)\cdot H(z) $ Discrete tustin', linewidth=1)
    ax1_impulse.plot(h_n_zoh[0], h_n_zoh[1][0], "C2", alpha=1, label=r'$ U(z)\cdot H(z)$ Discrete zoh', linewidth=1)

    #ax1_impulse.axvline(1 / ts * np.pi, color='red', ls=":")  # cutoff frequency
    ax1_impulse.legend()
    #ax1_impulse.set_xlim(0.01, 400)
    #ax1_impulse.set_ylim(-80, 5)
    ax1_impulse.set_xlabel('Time [seconds]')
    ax1_impulse.set_ylabel('Amplitude')
    ax1_impulse.margins(0, 0.1)
    ax1_impulse.grid(which='both', axis='both')


    plot_step.tight_layout()


    return plot_bode, plot_step



# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    #print_hi('PyCharm')
    # Create the analog filters to compare
    lowpass = sig.butter(1, 1, "low", analog=True)
    highpass = sig.butter(1, 80, "high", analog=True)
    bandstop = sig.butter(1, [90, 110], "bandstop", analog=True)
    elliptic = sig.ellip(8, 3, 20, 10, analog=True)

    # Make the figures for the comparison between the analog and discrete systems
    bode_lowpass, step_lowpass = compare2disc(lowpass, ts=0.01, title="lowpass-filter",
                                              fig_num=[1, 2], step_time=[0, 4])
    bode_highpass, step_highpass = compare2disc(highpass, ts=0.01, title="highpass-filter",
                                                fig_num=[3, 4], step_time=[0, 1])
    bode_bandstop, step_bandstop = compare2disc(bandstop, ts=0.01, title="bandstop-filter",
                                                fig_num=[5, 6], step_time=[0, 2])
    bode_elliptic, step_elliptic = compare2disc(elliptic, ts=0.01, title="elliptic-filter",
                                                fig_num=[7, 8], step_time=[0, 10])

    # Display the figures if enabled
    if show_figures:
        bode_lowpass.show()
        step_lowpass.show()
        bode_highpass.show()
        step_highpass.show()
        bode_bandstop.show()
        step_bandstop.show()
        bode_elliptic.show()
        step_elliptic.show()

    # Save the figures if enabled
    if save_figures:
        bode_lowpass.savefig("bode_lowpass")
        step_lowpass.savefig("step_lowpass")
        bode_highpass.savefig("bode_highpass")
        step_highpass.savefig("step_highpass")
        bode_bandstop.savefig("bode_bandstop")
        step_bandstop.savefig("step_bandstop")
        bode_elliptic.savefig("bode_elliptic")
        step_elliptic.savefig("step_elliptic")


# See PyCharm help at https://www.jetbrains.com/help/pycharm/

