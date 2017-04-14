# Trusted computing

Trusted computing is a technology that aims to enforce that computers consistently behave in expected ways. It aims to make computers safer and less prone to viruses and malware. This is achieved by using a hardware chip called Trusted Platform Module (TPM). Each TPM has a unique encryption key which cannot be accessed from the rest of the system.

# Trusted Platform Module (TPM)

TPM is a microchip on the motherboard of most contemporary PCs and laptops. It stores cryptographic keys and other sensitive data in its shielded memory and provides waysto securely use these keys. It supports a small set of cryptographic functions like key generation, signing, encryption, hashing, and MAC.

The core of TPM are the following three functionalities:

 * secure storage where the user can store data that is encrypted by keys available only to the TPM
 * platform integrity measurement and reporting which enables the platform to create reports of its integrity and configuration state 
 * platform authentication

TPM must be trusted by all parties to behave correctly and not being compromised - it acts as a Root of Trust.

## TPM initialization

Initialization of TPM requires physical access to the computer to turn on the TPM (to prevent malicious software being able to initialize TPM) - BIOS screen prompts and required strokes may be different on different computers.

In a TPM initialization wizard an owner password needs to be created (taking ownership of a TPM). The owner of the TPM can make a full use of TMP capabilities. Only the TPM owner can enable, disable, or clear the TPM without pyhsical access to the computer.

## HMAC in TPM

A TPM object (for example signing key) can have an associated authorization value that is known to the TPM and an authorized user of the key. The user prepares a TPM command (which includes the above-mentioned TPM object) and calculates an HMAC over the command using an HMAC key derived from the authorization value.

Sometimes HMAC key may be known only to TPM as well - for example when TPM outputs a value to be stored externally for later use, it includes an HMAC which can be later verified (this is called a ticket).

# Verified boot

Flash memory is a non-volatile (it is not affected by turning off and on again) storage that can be electrically erased and reprogrammed. It is used in memory cards, USB flash drives, solid-state drives (it uses transistors with two gates - electricity stays between the gates when the power is turned off).

Firmware is held in non-volatile memory devices such as ROM, EPROM, and flash memory. Firmware might never change during its economic lifetime. Changing/updating firmware may require physical replacement of integrated circuits or flash memory to be reprogrammed through a special process. 

Since the appearance of EEPROMs and flash memory, firmware usually does not reside in an unwriteable ROM, which means it can be changed (and what is making the distinction between firmware and software blurry). The distinction can be seen as: sofware is stored on a disk, firmware on flash memory.

BIOS is a type of firmware. It the first software that is run when a computer is powered on. BIOS initializes the hardware components and provides an abstraction layer for the hardware - consistent way for operation system and applications to interact with the keyboard, dispaly and other IO devices. Originally, BIOS was stored in a ROM chip, nowadays it is on flash memory so it can be updated without removing the chip from the motherboard.

BIOS is the first code that runs on CPU. If there is no integrity checking, BIOS is a great place to insert backdoors.

## Chromebook verified boot

Two pieces of hardware are used for security in Chromebook: a custom firmware chip and TPM. The customer firmware chip consists of read-only firmware and a read-write firmware that can be updated.

When the Chromebook is started (power on), a read-only firmware starts a verified boot - it reads the read-write firmware (let's denote it as blob F) and verifies it (Google generated a signature over hash(F) which is verified by using a public key stored in read-only memory).

If signature is correct, read-write firmware is run which performs a similar verification operation on the operating system kernel. The operating system kernel then continue the verification process for the system software. 
This way it is cryptographically assured that system code has not been modified by an attacker.
To prevent rollback attacks outdated signatures are stored in NVRAM in TPM to ensure that these are not accepted. Note that NVRAM location is lockable, such that only the read-only firmware can write to it.

# ARM TrustZone

...

# Intel TXT

## CPU internals

CPU is split into four different rings to provide separation of the OS and the applications that are managed by the OS.

It seems that it was intended that ring 0 is for OS, ring 1 for drivers, ring 2 for services, and ring 3 for applications. However, OS vendors used only ring 0 and ring 3 (all system related activities are in ring 0, instead of rings 0-2).

Allowing drivers to access ring 0 leads into many attacks.

CPU has a number of registers (for example 32) which are for storing the numbers CPU calculates with. CPU registers are a small amount of data storage which are directly accessible from CPU - no other memory is faster.
Electrical signals travel near light speed which is fast, but still finite (it could be something like 7 cm per clock tick), thus the memory needs to be as close as possible to the execution engine.

Every processor has also L1 cache and L2 cache, located further from the core. Even further away is RAM.

## Memory manipulation

Ring structure provides a protection to keep programs from accessing each other's memory. However, programs that run at ring 0 have access to all memory (including ring 3).
The easiest way to get an access to ring 0, is to install a driver (compromising OS is possible, but harder).
The only defense is to ensure that nothing happens to the code in ring 0.

Memory is also accessable using Direct Memory Access (DMA). DMA allows devices (like graphic, network, and sound cards) to access memory (all of it) without CPU involvement (thus it avoids all protections provided by CPU).
With DMA, the CPU can execute other operations while waiting for the DMA interrupt (notifying CPU that the transfer is done).

The attacker can convince the device to perform a DMA to the memory he is interested in.

## What is the purpose of TXT

The purpose of Intel TXT is to provide a trusted way for loading and executing system software (like OS kernel or Virtualization Machine Monitor (VMM)). To achieve this software measurements are performed and stored in TPM registers. 

Intel TXT is supposed to load OS kernel or VMM in secure way even if the boot sector is compromised.

TXT load process is called Late Launch and is implemented via a special CPU instruction named SENTER (AMD has its own version of late launch which is implemented using SKINIT instruction).

TXT provides only launch-time protection (not runtime) - it ensures that the code we load is what we really intended to load.

## Attackers considered by TXT

**TXT does not try to prevent the malicious software to be installed.** The attacker has the ability to run software processes on the platform. Hardware attacks are not considered (hardware attacks require physical access and are thus less likely).

When data is under TXT protection, no program (in whatever ring) can be able to view or modify it. This does not mean only protecting a physical memory page, but also resources that are attached to the application (that handles data): CPU registers, threads, debug operations, CPU counters (average rate per second at which context switches among threads, average rate per second at which the processor handles interrupts, instantaneous count of threads that are in the processor queue...).

## TXT features

Protected execution provides isolation of all resources associated with an application.

Protected memory pages - to protect against DMA from outside devices, MCH has the NoDMA table which specifies which memory pages are under protection of the CPU. When any device other than CPU operating in protected execution attempts to access a memory page, the access is blocked.

Sealed storage - data held by one configuration is not available to any other configuration.

Protected input - TXT provides a mechanism to establish a trusted channel between the protected execution environment and the keyboard and mouse.

Protected graphics - trusted channel between the application and the display adapter.

Attestations - to guarantee that a remote system is using TXT and for reporting how the protected execution environment is currently executing.


