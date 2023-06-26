#Importing Dependencies
from plotting import plot_xvg, CMD
import shlex, os


def run_cmd(command: CMD):
    match command:
        case CMD(command="help" | "h"):
            print('This is a help bar.')

        case CMD(command="gmplot", arguments=["plot","-f" | "--file" , file, *rest]):
            try:
                print(rest)
                plot_xvg(os.path.abspath(file), others=rest).runplot()

            except FileNotFoundError:
                print(f"'{file}' does not exist.")

        case CMD(command= "quit" | "exit" | "bye" | "end" | "q"):
            quit()

        case _:
            print(f"Unknown command: {str(command)!r}.")


def main() -> None:
    while True:
        cmd, *arguments = shlex.split(input("gmx+plotter >: $ "))
        run_cmd(CMD(cmd, arguments))


if __name__ == "__main__":
    main()
