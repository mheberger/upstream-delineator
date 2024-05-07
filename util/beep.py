import winsound


# Your script code here
def beep():
    # Play a beep sound when done
    frequency = 600
    duration = 1000
    winsound.Beep(frequency, duration)


if __name__ == "__main__":
    beep()